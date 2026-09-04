"""
plot_trace.py
--------------
Sanger kromatogramını çizer ve "mixed_candidate" olarak işaretlenmiş
pozisyonları dikey çizgilerle vurgular. Ayrıca her şüpheli pozisyon için
yakınlaştırılmış (zoom) bir alt-grafik üretir, böylece biyolog sadece o
birkaç pozisyona bakar.
"""

import matplotlib
matplotlib.use("Agg")  # dosyaya kaydetmek için, ekran gerektirmez
import matplotlib.pyplot as plt

from .abif_reader import SangerTrace

CHANNEL_COLORS = {"A": "#2E7D32", "C": "#1565C0", "G": "#F9A825", "T": "#C62828"}


def plot_full_trace(trace: SangerTrace, df, out_path: str, title: str = ""):
    fig, ax = plt.subplots(figsize=(16, 4))
    for base, color in CHANNEL_COLORS.items():
        ax.plot(trace.channels.get(base, []), color=color, linewidth=0.7, label=base)

    flagged = df[df["flag"] == "mixed_candidate"]
    for _, row in flagged.iterrows():
        ax.axvline(row["scan_pos"], color="black", alpha=0.25, linewidth=1)

    ax.set_xlabel("Scan pozisyonu")
    ax.set_ylabel("Sinyal yoğunluğu")
    ax.set_title(title or f"Sanger kromatogramı — {len(flagged)} şüpheli pozisyon işaretli")
    ax.legend(loc="upper right", ncol=4, fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_flagged_zooms(trace: SangerTrace, df, out_path: str, pad: int = 20):
    flagged = df[df["flag"] == "mixed_candidate"]
    if flagged.empty:
        return None

    n = len(flagged)
    ncols = min(4, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3 * nrows), squeeze=False)

    for ax_idx, (_, row) in enumerate(flagged.iterrows()):
        ax = axes[ax_idx // ncols][ax_idx % ncols]
        center = row["scan_pos"]
        lo, hi = max(0, center - pad), center + pad
        for base, color in CHANNEL_COLORS.items():
            channel = trace.channels.get(base, [])
            ax.plot(range(lo, min(hi, len(channel))), channel[lo:hi], color=color, linewidth=1.2)
        ax.axvline(center, color="black", alpha=0.3, linestyle="--")
        ax.set_title(
            f"poz {row['index']} ({row['called_base']}): "
            f"{row['primary_base']}/{row['secondary_base']} oran={row['ratio']}",
            fontsize=9,
        )
        ax.set_xticks([])

    # kullanılmayan alt grafikleri gizle
    for ax_idx in range(n, nrows * ncols):
        axes[ax_idx // ncols][ax_idx % ncols].axis("off")

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    return out_path
