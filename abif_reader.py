"""
abif_reader.py
----------------
Applied Biosystems ABIF (.ab1) dosya formatı için bağımsız (Biopython'sız) okuyucu.

ABIF formatı, ABI genetik analiz cihazlarının (örn. Sanger dizileme cihazları)
ürettiği ikili (binary) dosya formatıdır. Format kamuya açık şekilde
dokümante edilmiştir (Applied Biosystems ABIF File Format spec).

Yapı:
  - Dosya başında 4 byte 'ABIF' + 2 byte versiyon
  - Ardından tek bir "root" dizin girişi (28 byte) -> bu girişin
    dataoffset alanı, gerçek dizin girişlerinin (her biri 28 byte)
    bulunduğu konumu gösterir.
  - Her dizin girişi bir "tag" (örn. DATA9, PLOC1, PBAS1, FWO_1) tanımlar
    ve bu tag'e ait verinin dosyadaki konumunu / tipini / eleman sayısını verir.

Bu modül sadece bizim ihtiyacımız olan tag'leri (DATA9-12, PLOC1/2,
PBAS1/2, PCON1/2, FWO_1) okur ama genel bir sözlük (dict) olarak
tüm tag'leri de sunar.
"""

import struct
from dataclasses import dataclass, field


DIR_ENTRY_FMT = ">4siHHiiii"  # name(4s) number(i) elementtype(h) elementsize(h) numelements(i) datasize(i) dataoffset(i) datahandle(i)
DIR_ENTRY_SIZE = struct.calcsize(DIR_ENTRY_FMT)  # 28 byte olmalı


class AbifFormatError(Exception):
    """Dosya ABIF formatına uymuyorsa fırlatılır."""


def _parse_value(raw: bytes, elementtype: int, numelements: int):
    """Ham byte dizisini elementtype'a göre Python nesnesine çevirir."""
    if elementtype == 2:  # char array (ASCII metin)
        return raw[:numelements].decode("ascii", errors="replace")
    if elementtype == 3:  # word (uint16)
        return list(struct.unpack(f">{numelements}H", raw[: 2 * numelements]))
    if elementtype == 4:  # short (int16)
        return list(struct.unpack(f">{numelements}h", raw[: 2 * numelements]))
    if elementtype == 5:  # long (int32)
        return list(struct.unpack(f">{numelements}i", raw[: 4 * numelements]))
    if elementtype == 7:  # float
        return list(struct.unpack(f">{numelements}f", raw[: 4 * numelements]))
    if elementtype == 8:  # double
        return list(struct.unpack(f">{numelements}d", raw[: 8 * numelements]))
    if elementtype == 10:  # date (year:int16, month:byte, day:byte)
        return raw[:numelements]
    if elementtype == 11:  # time
        return raw[:numelements]
    if elementtype == 13:  # bool
        return bool(raw[0]) if raw else None
    if elementtype == 18:  # pascal string: ilk byte uzunluk
        length = raw[0]
        return raw[1 : 1 + length].decode("ascii", errors="replace")
    if elementtype == 19:  # C string (null-terminated)
        end = raw.find(b"\x00")
        return raw[: end if end != -1 else len(raw)].decode("ascii", errors="replace")
    # Bilinmeyen / kullanıcı tanımlı tipler: ham byte olarak döndür
    return raw


def read_abif(path: str) -> dict:
    """
    .ab1 dosyasını okur ve {tag_adi+tag_no: değer} şeklinde bir dict döndürür.
    Örnek anahtarlar: 'DATA9', 'DATA10', 'DATA11', 'DATA12', 'PLOC1', 'PLOC2',
                       'PBAS1', 'PBAS2', 'PCON1', 'PCON2', 'FWO_1'
    """
    with open(path, "rb") as f:
        data = f.read()

    if len(data) < 30 or data[0:4] != b"ABIF":
        raise AbifFormatError(
            f"'{path}' geçerli bir ABIF (.ab1) dosyası gibi görünmüyor "
            f"(dosya başında 'ABIF' imzası bulunamadı)."
        )

    # Root dizin girişi, offset 6'dan başlar (4 byte magic + 2 byte version)
    root_raw = data[6 : 6 + DIR_ENTRY_SIZE]
    (_name, _number, _etype, _esize, numelements, _datasize, dataoffset, _dhandle) = (
        struct.unpack(DIR_ENTRY_FMT, root_raw)
    )

    tags = {}
    for i in range(numelements):
        entry_offset = dataoffset + i * DIR_ENTRY_SIZE
        entry_raw = data[entry_offset : entry_offset + DIR_ENTRY_SIZE]
        if len(entry_raw) < DIR_ENTRY_SIZE:
            break
        name, number, etype, esize, nelem, dsize, doffset, dhandle = struct.unpack(
            DIR_ENTRY_FMT, entry_raw
        )
        name = name.decode("ascii", errors="replace").strip()

        # datasize <= 4 byte ise veri, dataoffset alanının kendisinde (inline) saklanır
        if dsize <= 4:
            raw = entry_raw[20:24]
        else:
            raw = data[doffset : doffset + dsize]

        try:
            value = _parse_value(raw, etype, nelem)
        except struct.error:
            value = raw  # parse edilemeyen nadir tag'leri ham bırak

        tags[f"{name}{number}"] = value

    return tags


@dataclass
class SangerTrace:
    """Bir .ab1 dosyasından çıkarılan, analiz için gereken temel veriler."""

    channel_order: str  # örn. "GATC" -> DATA9,10,11,12 hangi baza karşılık geliyor
    channels: dict  # {'A': [int,...], 'C': [...], 'G': [...], 'T': [...]}
    bases: str  # basecaller'ın okuduğu dizi (örn. "ACGGT...")
    peak_locations: list  # her baz için trace üzerindeki scan pozisyonu
    qualities: list = field(default_factory=list)  # her baz için kalite skoru (varsa)
    sample_name: str = ""
    source_file: str = ""


def load_sanger_trace(path: str) -> SangerTrace:
    """
    .ab1 dosyasını okuyup analiz için hazır bir SangerTrace nesnesi döndürür.
    Basecaller'ın orijinal (edit edilmemiş) çağrılarını kullanır (tag no=1);
    çoğu ham .ab1 dosyasında sadece bunlar mevcuttur.
    """
    tags = read_abif(path)

    channel_order = tags.get("FWO_1", "GATC")
    if not isinstance(channel_order, str) or len(channel_order) != 4:
        channel_order = "GATC"  # ABI cihazlarının fabrika varsayılanı

    data_tags = ["DATA9", "DATA10", "DATA11", "DATA12"]
    missing = [t for t in data_tags if t not in tags]
    if missing:
        raise AbifFormatError(
            f"Beklenen trace kanalları ({', '.join(missing)}) dosyada bulunamadı. "
            "Bu dosya standart bir Sanger .ab1 dosyası olmayabilir."
        )

    channels = {base: tags[tag] for base, tag in zip(channel_order, data_tags)}

    # Edit edilmiş (2) sürüm varsa onu tercih et, yoksa orijinal (1) sürümü kullan
    bases = tags.get("PBAS2") or tags.get("PBAS1") or ""
    peak_locations = tags.get("PLOC2") or tags.get("PLOC1") or []
    qualities = tags.get("PCON2") or tags.get("PCON1") or []

    sample_name = tags.get("SMPL1", "") if isinstance(tags.get("SMPL1", ""), str) else ""

    return SangerTrace(
        channel_order=channel_order,
        channels=channels,
        bases=bases,
        peak_locations=list(peak_locations),
        qualities=list(qualities),
        sample_name=sample_name,
        source_file=path,
    )
