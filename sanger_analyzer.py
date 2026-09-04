"""
sanger_analyzer.py
--------------------
Sanger dizileme (.ab1) trace verisinde her baz pozisyonu için 4 kanalın
(A, C, G, T) pik yüksekliklerini karşılaştırıp, "tek pik" (homozigot/normal)
ile "çift pik" (heterozigot/mixed/mutasyon şüphesi) pozisyonlarını ayırt eder.

Mantık (klinik sekans analiz yazılımlarında -Mutation Surveyor, 4Peaks vb.-
kullanılan standart yaklaşıma benzer):
  1. Basecaller'ın işaretlediği her pozisyonun (PLOC) etrafında küçük bir
     pencere içinde her kanalın maksimum yüksekliğine bakılır (basecaller'ın
     bulduğu tam pikin merkezi birkaç sample kayabilir).
  2. 4 kanal yüksekliğe göre sıralanır: en yüksek = "primary", ikinci en
     yüksek = "secondary".
  3. secondary / primary oranı (peak ratio) bir eşik değerin (varsayılan
     %25) üzerindeyse ve secondary sinyali arka plan gürültüsünün belirgin
     şekilde üzerindeyse -> pozisyon "MIXED / heterozigot şüphesi" olarak
     işaretlenir.
  4. İsteğe bağlı referans (yabani tip / wild-type) dizisi verilirse, basit
     hizalama ile hangi pozisyonların referanstan farklı olduğu da eklenir.

Bu eşik değerleri laboratuvar/panel'e göre ayarlanabilir olacak şekilde
parametrik bırakılmıştır -- gerçek örnek verileriyle kalibre edilmesi önerilir.
"""

from dataclasses import dataclass
from typing import Optional
import statistics

import pandas as pd

from .abif_reader import SangerTrace, load_sanger_trace

BASES = ("A", "C", "G", "T")


@dataclass
class PeakCall:
    index: int  # dizideki pozisyon (0-tabanlı)
    scan_pos: int  # trace üzerindeki ham scan koordinatı
    called_base: str  # basecaller'ın orijinal çağrısı
    primary_base: str
    primary_height: int
    secondary_base: str
    secondary_height: int
    ratio: float  # secondary / primary
    flag: str  # "normal" | "mixed_candidate" | "low_signal"
    ref_base: Optional[str] = None
    matches_ref: Optional[bool] = None


def _local_max(channel: list, center: int, window: int) -> int:
    """center etrafında +/-window örneklem içindeki maksimum yüksekliği döndürür."""
    lo = max(0, center - window)
    hi = min(len(channel), center + window + 1)
    if lo >= hi:
        return 0
    return max(channel[lo:hi])


def analyze_trace(
    trace: SangerTrace,
    window: int = 4,
    mixed_ratio_threshold: float = 0.25,
    noise_floor_k: float = 3.0,
    reference: Optional[str] = None,
) -> pd.DataFrame:
    """
    Her baz pozisyonu için pik analizini yapar ve bir pandas DataFrame döndürür.

    Parametreler
    ----------
    window : her PLOC pozisyonu etrafında pik aranacak sample penceresi
    mixed_ratio_threshold : secondary/primary oranı bu değerin üzerindeyse
        pozisyon "mixed_candidate" olarak işaretlenir (varsayılan 0.25 ~
        literatürde heterozigot tespiti için sık kullanılan eşik).
    noise_floor_k : arka plan gürültüsünü (std) kaç katı aşan sinyaller
        "gerçek" sayılır; bunun altındaki secondary pikler gürültü sayılıp
        işaretlenmez.
    reference : varsa, karşılaştırma için referans (yabani tip) dizisi.
        Basit indeks-bazlı hizalama kullanılır (insersiyon/delesyon
        olmayan, tek nokta mutasyonu taramaları için uygundur).
    """
    if not trace.bases or not trace.peak_locations:
        raise ValueError(
            "Trace içinde baz çağrısı (PBAS) veya pik konumu (PLOC) bulunamadı."
        )

    # Genel gürültü tabanını tahmin etmek için medyan + MAD (median absolute
    # deviation) kullanılır. Trace'in büyük çoğunluğu piklerin ARASINDAKİ
    # düşük seviyeli gürültüden oluştuğu için medyan/MAD, ortalama/std'ye göre
    # çok daha sağlamdır (ortalama, yüksek piklerden kolayca şişer ve gerçek
    # ikincil pikleri "gürültü" gibi gösterip kaçırmaya sebep olur).
    all_signal = [v for ch in trace.channels.values() for v in ch]
    if all_signal:
        med = statistics.median(all_signal)
        mad = statistics.median(abs(v - med) for v in all_signal) or 1.0
    else:
        med, mad = 0, 1.0
    noise_floor = med + noise_floor_k * mad

    rows = []
    n = min(len(trace.bases), len(trace.peak_locations))
    for i in range(n):
        called_base = trace.bases[i].upper()
        scan_pos = trace.peak_locations[i]

        heights = {}
        for base in BASES:
            channel = trace.channels.get(base, [])
            heights[base] = _local_max(channel, scan_pos, window)

        ranked = sorted(heights.items(), key=lambda kv: kv[1], reverse=True)
        (primary_base, primary_height), (secondary_base, secondary_height) = ranked[0], ranked[1]

        ratio = (secondary_height / primary_height) if primary_height > 0 else 0.0

        if primary_height < noise_floor:
            flag = "low_signal"  # bu pozisyondaki genel sinyal çok zayıf, dikkatli yorumla
        elif ratio >= mixed_ratio_threshold and secondary_height > noise_floor:
            flag = "mixed_candidate"  # olası heterozigot / iki pik / mutasyon şüphesi
        else:
            flag = "normal"

        row = PeakCall(
            index=i,
            scan_pos=scan_pos,
            called_base=called_base,
            primary_base=primary_base,
            primary_height=primary_height,
            secondary_base=secondary_base,
            secondary_height=secondary_height,
            ratio=round(ratio, 3),
            flag=flag,
        )

        if reference and i < len(reference):
            row.ref_base = reference[i].upper()
            row.matches_ref = row.ref_base == primary_base

        rows.append(row)

    df = pd.DataFrame([r.__dict__ for r in rows])
    return df


def summarize(df: pd.DataFrame) -> dict:
    """Hızlı bir özet: kaç pozisyon var, kaçı flagli, vb."""
    counts = df["flag"].value_counts().to_dict()
    summary = {
        "toplam_pozisyon": len(df),
        "normal": counts.get("normal", 0),
        "mixed_candidate": counts.get("mixed_candidate", 0),
        "low_signal": counts.get("low_signal", 0),
    }
    if "matches_ref" in df.columns and df["matches_ref"].notna().any():
        summary["referanstan_farkli"] = int((df["matches_ref"] == False).sum())  # noqa: E712
    return summary


def flagged_positions(df: pd.DataFrame) -> pd.DataFrame:
    """Sadece dikkat gerektiren (mixed_candidate) pozisyonları döndürür -
    biyoloğun her pike değil, sadece bu satırlara bakması yeterli olur."""
    return df[df["flag"] == "mixed_candidate"].copy()


def analyze_file(
    path: str,
    window: int = 4,
    mixed_ratio_threshold: float = 0.25,
    noise_floor_k: float = 3.0,
    reference: Optional[str] = None,
) -> pd.DataFrame:
    """Tek adımda: dosyayı oku + analiz et."""
    trace = load_sanger_trace(path)
    return analyze_trace(
        trace,
        window=window,
        mixed_ratio_threshold=mixed_ratio_threshold,
        noise_floor_k=noise_floor_k,
        reference=reference,
    )
