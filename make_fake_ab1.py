"""
make_fake_ab1.py
------------------
Gerçek bir Sanger cihazı çıktısı olmadan abif_reader.py ve sanger_analyzer.py'ı
uçtan uca test edebilmek için, bilinen bir dizi + bilinçli olarak yerleştirilmiş
"çift pik" (heterozigot) pozisyonları içeren SENTETİK bir .ab1 dosyası üretir.

Bu dosya gerçek biyolojik veri DEĞİLDİR, sadece format/algoritma doğrulaması içindir.
"""

import struct
import random

DIR_ENTRY_FMT = ">4siHHiiii"
DIR_ENTRY_SIZE = struct.calcsize(DIR_ENTRY_FMT)


def make_gaussian_trace(length, peak_positions, peak_heights, width=3):
    """Basit gaussian-benzeri pikler içeren sahte bir trace üretir."""
    trace = [0.0] * length
    for pos, height in zip(peak_positions, peak_heights):
        for x in range(max(0, pos - 4 * width), min(length, pos + 4 * width)):
            trace[x] += height * pow(2.71828, -((x - pos) ** 2) / (2 * width * width))
    return [int(min(max(v, 0), 32000)) for v in trace]


def build_fake_ab1(out_path, sequence, mixed_positions, spacing=20, noise=15, seed=42):
    """
    sequence: temel (primary) baz dizisi, örn. "ACGGT..."
    mixed_positions: {index: (secondary_base, secondary_height_ratio)} -> bu
        pozisyonlarda ikinci bir pik eklenir (heterozigot simülasyonu)
    """
    random.seed(seed)
    n = len(sequence)
    peak_locations = [spacing * (i + 1) for i in range(n)]
    trace_length = peak_locations[-1] + spacing

    base_order = "GATC"
    channels = {b: [0] * n for b in "ACGT"}  # her pozisyon için yükseklik (sonra trace'e çevrilecek)

    primary_positions = {b: [] for b in "ACGT"}
    primary_heights = {b: [] for b in "ACGT"}

    for i, base in enumerate(sequence):
        base = base.upper()
        height = random.randint(800, 1200)
        primary_positions[base].append(peak_locations[i])
        primary_heights[base].append(height)

        if i in mixed_positions:
            sec_base, ratio = mixed_positions[i]
            sec_height = int(height * ratio)
            primary_positions[sec_base].append(peak_locations[i])
            primary_heights[sec_base].append(sec_height)

    channel_traces = {}
    for b in "ACGT":
        tr = make_gaussian_trace(trace_length, primary_positions[b], primary_heights[b])
        # arka plan gürültüsü ekle
        tr = [max(0, v + random.randint(0, noise)) for v in tr]
        channel_traces[b] = tr

    # DATA9-12 sırası FWO_ ile eşleşmeli: base_order[k] -> DATAtag[k]
    data_tag_numbers = [9, 10, 11, 12]
    data_arrays = [channel_traces[base_order[k]] for k in range(4)]

    tags_to_write = []  # (name, number, elementtype, numelements, raw_bytes)

    def short_array_tag(name, number, values):
        raw = struct.pack(f">{len(values)}h", *values)
        return (name, number, 4, len(values), raw)

    def char_tag(name, number, text):
        raw = text.encode("ascii")
        return (name, number, 2, len(raw), raw)

    for k in range(4):
        tags_to_write.append(short_array_tag("DATA", data_tag_numbers[k], data_arrays[k]))

    tags_to_write.append(short_array_tag("PLOC", 1, peak_locations))
    tags_to_write.append(char_tag("PBAS", 1, sequence))
    tags_to_write.append(short_array_tag("PCON", 1, [40] * n))
    tags_to_write.append(char_tag("FWO_", 1, base_order))
    tags_to_write.append(char_tag("SMPL", 1, "TEST_SAMPLE_01"))

    # --- Dosyayı bayt bayt inşa et ---
    header = b"ABIF" + struct.pack(">h", 101)  # magic + version
    root_entry_offset = 6
    dir_array_offset = root_entry_offset + DIR_ENTRY_SIZE  # dizin girişleri hemen sonra başlasın

    num_tags = len(tags_to_write)
    dir_entries_bytes = b""
    data_blob = b""
    data_blob_base_offset = dir_array_offset + num_tags * DIR_ENTRY_SIZE

    for name, number, etype, nelem, raw in tags_to_write:
        datasize = len(raw)
        if datasize <= 4:
            padded = raw + b"\x00" * (4 - datasize)
            entry = struct.pack(
                DIR_ENTRY_FMT,
                name.encode("ascii"), number, etype, 2 if etype == 2 else 2,
                nelem, datasize, struct.unpack(">i", padded)[0], 0,
            )
        else:
            offset = data_blob_base_offset + len(data_blob)
            entry = struct.pack(
                DIR_ENTRY_FMT,
                name.encode("ascii"), number, etype, 2 if etype == 2 else 2,
                nelem, datasize, offset, 0,
            )
            data_blob += raw
        dir_entries_bytes += entry

    root_entry = struct.pack(
        DIR_ENTRY_FMT,
        b"tdir", 1, 1023, DIR_ENTRY_SIZE, num_tags,
        num_tags * DIR_ENTRY_SIZE, dir_array_offset, 0,
    )

    full = header + root_entry + dir_entries_bytes + data_blob

    with open(out_path, "wb") as f:
        f.write(full)

    return out_path


if __name__ == "__main__":
    seq = "ACGGTTCAGCTAGGCATTAGCCA"
    # index 8 (base 'G') üzerine %40 oranında 'A' ikinci piki ekle -> heterozigot simülasyonu
    # index 15 (base 'G') üzerine %60 oranında 'T' ikinci piki ekle -> güçlü heterozigot
    mixed = {8: ("A", 0.40), 15: ("T", 0.60)}
    path = build_fake_ab1("/home/claude/sanger_tool/test_sample.ab1", seq, mixed)
    print(f"Sahte .ab1 dosyası oluşturuldu: {path}")
    print(f"Dizi: {seq}")
    print(f"Beklenen mixed pozisyonlar: {mixed}")
