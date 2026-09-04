"""
main.py
--------
Kullanım örnekleri:

  python3 main.py ornek.ab1
  python3 main.py ornek.ab1 --ref ACGTACGT... --out rapor.csv --plot kromatogram.png
  python3 main.py ornek.ab1 --threshold 0.20

Çıktı:
  - Terminale özet (kaç pozisyon normal / mixed_candidate / low_signal)
  - CSV raporu (her pozisyonun detayı; sadece flag='mixed_candidate' olanlar
    "onay gerektiren" pozisyonlardır)
  - (opsiyonel) Kromatogram görseli + şüpheli pozisyonların zoom görünümü
"""

import argparse
import sys

from .abif_reader import load_sanger_trace, AbifFormatError
from .sanger_analyzer import analyze_trace, summarize, flagged_positions
from .plot_trace import plot_full_trace, plot_flagged_zooms


def main():
    parser = argparse.ArgumentParser(description="Sanger (.ab1) kromatogram analiz aracı")
    parser.add_argument("ab1_path", help=".ab1 dosyasının yolu")
    parser.add_argument("--ref", default=None, help="Referans (yabani tip) dizisi (opsiyonel)")
    parser.add_argument("--out", default="rapor.csv", help="CSV rapor çıktı yolu")
    parser.add_argument("--plot", default=None, help="Kromatogram görseli için PNG yolu (opsiyonel)")
    parser.add_argument("--threshold", type=float, default=0.25, help="Mixed pik oranı eşiği (varsayılan 0.25)")
    parser.add_argument("--window", type=int, default=4, help="Pik arama penceresi, sample cinsinden")
    args = parser.parse_args()

    try:
        trace = load_sanger_trace(args.ab1_path)
    except AbifFormatError as e:
        print(f"HATA: {e}", file=sys.stderr)
        sys.exit(1)

    df = analyze_trace(
        trace,
        window=args.window,
        mixed_ratio_threshold=args.threshold,
        reference=args.ref,
    )

    df.to_csv(args.out, index=False)

    summary = summarize(df)
    print(f"\nDosya: {args.ab1_path}")
    print(f"Toplam pozisyon: {summary['toplam_pozisyon']}")
    print(f"  Normal (tek pik): {summary['normal']}")
    print(f"  ŞÜPHELİ (mixed/heterozigot adayı): {summary['mixed_candidate']}  <-- buraya bak")
    print(f"  Düşük sinyal: {summary['low_signal']}")
    if "referanstan_farkli" in summary:
        print(f"  Referanstan farklı: {summary['referanstan_farkli']}")
    print(f"\nDetaylı rapor kaydedildi: {args.out}")

    flagged = flagged_positions(df)
    if not flagged.empty:
        print("\nŞüpheli pozisyonlar:")
        print(flagged[["index", "called_base", "primary_base", "secondary_base", "ratio"]].to_string(index=False))

    if args.plot:
        plot_full_trace(trace, df, args.plot)
        print(f"\nKromatogram görseli kaydedildi: {args.plot}")
        zoom_path = args.plot.replace(".png", "_zoom.png")
        if plot_flagged_zooms(trace, df, zoom_path):
            print(f"Şüpheli pozisyonların yakın görünümü: {zoom_path}")


if __name__ == "__main__":
    main()
