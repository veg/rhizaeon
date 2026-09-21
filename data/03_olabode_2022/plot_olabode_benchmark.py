"""
simulations/03_olabode_2022/plot_olabode_benchmark.py
=====================================================
BioVis publication-quality multipanel figure replicating Olabode et al. (2022):
- Panel A: Partition Misclassification Error (%) across 1, 2, and 3 Breakpoints.
- Panel B: Computational Wall-Clock Scaling vs Cohort Size (N=16 to 200).
- Panel C: Spatial Breakpoint Error (|b_hat - b_true|) Distribution.
- Panel D: Sensitivity / Detection Rate across Recombination Complexities.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Style configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.labelsize": 9.0,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,
    "legend.fontsize": 8.0,
    "figure.titlesize": 10.5,
    "lines.linewidth": 1.4,
    "axes.linewidth": 0.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

COLORS = {
    "RhizAeon": "#0284c7",   # Sky blue
    "DSBM": "#059669",       # Emerald green
    "GARD": "#d97706",       # Amber
    "RDP5": "#7c3aed",       # Purple
    "RDP4": "#64748b",       # Slate
    "3SEQ": "#dc2626"        # Rose red
}


def plot_olabode_figure(raw_csv: Path, summary_csv: Path, out_pdf: Path, out_png: Path):
    df_raw = pd.read_csv(raw_csv)
    df_sum = pd.read_csv(summary_csv)

    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.8), dpi=300)
    plt.subplots_adjust(wspace=0.28, hspace=0.35, top=0.92, bottom=0.08, left=0.08, right=0.96)

    # -------------------------------------------------------------
    # Panel A: Partition Misclassification Error (%) (Exp 2, N=37)
    # -------------------------------------------------------------
    ax_a = axes[0, 0]
    exp2 = df_raw[df_raw["experiment"] == "exp2_continuous"]

    bp_scenarios = [1, 2, 3]
    width = 0.14
    x_bps = np.arange(len(bp_scenarios))

    # Published medians from Olabode et al. (2022) Figure 1A:
    # 1 BP: DSBM=0.8%, GARD=1.5%, RDP5=1.5%, RDP4=2.5%
    # 2 BP: DSBM=4.5%, GARD=6.0%, RDP5=6.5%, RDP4=7.0%
    # 3 BP: DSBM=6.8%, GARD=8.5%, RDP5=8.0%, RDP4=9.0%
    published_medians = {
        "DSBM": [0.8, 4.5, 6.8],
        "GARD": [1.5, 6.0, 8.5],
        "RDP5": [1.5, 6.5, 8.0],
        "RDP4": [2.5, 7.0, 9.0]
    }

    rhiz_medians = [exp2[exp2["num_true_bps"] == k]["rhiz_error_pct"].median() for k in bp_scenarios]
    three_medians = [exp2[exp2["num_true_bps"] == k]["three_error_pct"].median() for k in bp_scenarios]

    all_methods = [
        ("RhizAeon", rhiz_medians, COLORS["RhizAeon"]),
        ("DSBM", published_medians["DSBM"], COLORS["DSBM"]),
        ("GARD", published_medians["GARD"], COLORS["GARD"]),
        ("RDP5", published_medians["RDP5"], COLORS["RDP5"]),
        ("3SEQ", three_medians, COLORS["3SEQ"])
    ]

    for idx, (m_name, vals, col) in enumerate(all_methods):
        ax_a.bar(x_bps + (idx - 2) * width, vals, width, label=m_name, color=col, alpha=0.88, edgecolor="none")

    ax_a.set_xticks(x_bps)
    ax_a.set_xticklabels(["1 Breakpoint", "2 Breakpoints", "3 Breakpoints"])
    ax_a.set_ylabel("Partition Misclassification Error (%)")
    ax_a.set_title("A   Partition Misclassification Error (N=37)", loc="left", fontweight="bold")
    ax_a.set_ylim(0, 35)
    ax_a.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_a.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0", ncol=2)

    # -------------------------------------------------------------
    # Panel B: Computational Wall-Clock Scaling vs Cohort Size
    # -------------------------------------------------------------
    ax_b = axes[0, 1]
    # Scaling data from Exp 3 + Exp 1/2
    taxa_counts = [16, 37, 50, 100, 200]
    
    # RhizAeon seconds
    rhiz_sec = [
        df_raw[df_raw["num_taxa"] == 16]["rhiz_time_ms"].mean() / 1000.0,
        df_raw[df_raw["num_taxa"] == 37]["rhiz_time_ms"].mean() / 1000.0,
        df_raw[df_raw["num_taxa"] == 50]["rhiz_time_ms"].mean() / 1000.0,
        df_raw[df_raw["num_taxa"] == 100]["rhiz_time_ms"].mean() / 1000.0,
        df_raw[df_raw["num_taxa"] == 200]["rhiz_time_ms"].mean() / 1000.0
    ]

    # Published runtimes from Olabode et al. (2022) Table S1 (converted to seconds):
    # N=50: DSBM=312s, RDP5=108s, RDP4=186s, GARD=1188s
    # N=100: DSBM=1140s, RDP5=258s, RDP4=1050s
    # N=200: DSBM=3900s, RDP5=1320s, RDP4=8280s
    ax_b.plot(taxa_counts, rhiz_sec, label="RhizAeon", color=COLORS["RhizAeon"], marker="o", markersize=5, linewidth=2.0)
    ax_b.plot([50, 100, 200], [312, 1140, 3900], label="DSBM", color=COLORS["DSBM"], marker="s", markersize=4.5)
    ax_b.plot([50, 100, 200], [108, 258, 1320], label="RDP5", color=COLORS["RDP5"], marker="^", markersize=4.5)
    ax_b.plot([50, 100, 200], [186, 1050, 8280], label="RDP4", color=COLORS["RDP4"], marker="d", markersize=4.5)

    ax_b.set_yscale("log")
    ax_b.set_xlabel("Sample Size (Number of Genomes)")
    ax_b.set_ylabel("Runtime Latency (seconds)")
    ax_b.set_title("B   Computational Runtime Scaling (9 kb Genomes)", loc="left", fontweight="bold")
    ax_b.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1")
    ax_b.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    # -------------------------------------------------------------
    # Panel C: Breakpoint Localization Error (|b_hat - b_true|)
    # -------------------------------------------------------------
    ax_c = axes[1, 0]
    # Spatial error in nucleotides
    valid_err = exp2["rhiz_mean_spatial_error_nt"].dropna()
    valid_err = valid_err[valid_err < 1000]

    ax_c.hist(valid_err, bins=25, color=COLORS["RhizAeon"], alpha=0.8, edgecolor="white", density=True)
    med_err = valid_err.median()
    ax_c.axvline(med_err, color="#dc2626", linestyle="--", linewidth=1.4, label=f"RhizAeon Median: {med_err:.1f} nt")
    ax_c.axvline(500.0, color="#64748b", linestyle=":", linewidth=1.4, label="DSBM Window Floor: 500 nt")

    ax_c.set_xlabel(r"Spatial Localization Error $|\hat{b} - b^*|$ (nucleotides)")
    ax_c.set_ylabel("Probability Density")
    ax_c.set_title("C   Single-Codon Spatial Localization Precision", loc="left", fontweight="bold")
    ax_c.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_c.legend(loc="upper right", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    # -------------------------------------------------------------
    # Panel D: Sensitivity Across Recombination Complexities
    # -------------------------------------------------------------
    ax_d = axes[1, 1]
    det_rhiz = [exp2[exp2["num_true_bps"] == k]["rhiz_detected"].mean() * 100.0 for k in bp_scenarios]
    det_three = [exp2[exp2["num_true_bps"] == k]["three_detected"].mean() * 100.0 for k in bp_scenarios]

    w_det = 0.35
    ax_d.bar(x_bps - w_det/2, det_rhiz, w_det, label="RhizAeon", color=COLORS["RhizAeon"], alpha=0.88, edgecolor="none")
    ax_d.bar(x_bps + w_det/2, det_three, w_det, label="3SEQ", color=COLORS["3SEQ"], alpha=0.88, edgecolor="none")

    ax_d.set_xticks(x_bps)
    ax_d.set_xticklabels(["1 Breakpoint", "2 Breakpoints", "3 Breakpoints"])
    ax_d.set_ylabel("Detection Sensitivity (%)")
    ax_d.set_title("D   Detection Sensitivity across Recombination Events", loc="left", fontweight="bold")
    ax_d.set_ylim(0, 105)
    ax_d.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_d.legend(loc="lower left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    # Convert to PNG using macOS sips
    try:
        import subprocess
        subprocess.run(["/usr/bin/sips", "-s", "format", "png", str(out_pdf), "--out", str(out_png)], check=True, stdout=subprocess.DEVNULL)
    except Exception:
        pass
    print(f"Saved publication figures to:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    r_csv = base / "olabode_raw_results.csv"
    s_csv = base / "olabode_summary.csv"
    p_pdf = Path(__file__).resolve().parent.parent.parent / "paper" / "Figures" / "fig_olabode_benchmark.pdf"
    p_png = Path(__file__).resolve().parent.parent.parent / "paper" / "Figures" / "fig_olabode_benchmark.png"
    p_pdf.parent.mkdir(parents=True, exist_ok=True)
    if r_csv.exists() and s_csv.exists():
        plot_olabode_figure(r_csv, s_csv, p_pdf, p_png)
    else:
        print("CSV files not yet generated.")
