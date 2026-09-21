"""
benchmarks/posada_2001_benchmark/plot_posada_benchmark.py
=========================================================
Generates publication-quality figure replicating Posada & Crandall (2001, PNAS Fig 1):
Left column: Power vs Recombination parameter rho across theta in {10, 50, 100, 200}.
Right column: False Positive Rate vs Rate Variation alpha across theta in {10, 50, 100, 200}.

Compares RhizAeon and 3SEQ directly against the historical benchmarks from Posada 2001.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Scientific styling
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 8.5,
    "figure.titlesize": 13,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.3,
    "lines.linewidth": 1.8,
    "lines.markersize": 5
})

# Historical benchmark values digitized from Posada & Crandall (2001) Figure 1
# for key reference methods: GENECONV, MAXCHI, CHIMAERA, HOMOPLASY TEST, RETICULATE, RDP, SIMPLOT
HISTORICAL_POWER = {
    # theta=10: rho in [0, 1, 4, 16, 64]
    10: {
        "GENECONV": [0, 5, 30, 48, 38],
        "MAXCHI": [0, 5, 22, 28, 28],
        "CHIMAERA": [0, 5, 32, 46, 38],
        "HOMOPLASY": [0, 15, 38, 82, 95],
        "RETICULATE": [0, 8, 25, 55, 72],
        "RDP": [0, 0, 2, 5, 8],
        "SIMPLOT": [0, 0, 0, 10, 20],
    },
    # theta=50: rho in [0, 1, 4, 16, 64]
    50: {
        "GENECONV": [0, 40, 78, 95, 96],
        "MAXCHI": [0, 30, 60, 92, 95],
        "CHIMAERA": [0, 40, 78, 95, 96],
        "HOMOPLASY": [0, 0, 12, 68, 88],
        "RETICULATE": [0, 28, 60, 92, 95],
        "RDP": [0, 15, 45, 82, 94],
        "SIMPLOT": [0, 0, 15, 40, 85],
    },
    # theta=100: rho in [0, 1, 4, 16, 64]
    100: {
        "GENECONV": [0, 45, 78, 96, 97],
        "MAXCHI": [0, 40, 70, 95, 96],
        "CHIMAERA": [0, 45, 75, 95, 97],
        "HOMOPLASY": [0, 0, 0, 35, 92],
        "RETICULATE": [0, 30, 65, 92, 95],
        "RDP": [0, 15, 55, 92, 96],
        "SIMPLOT": [0, 2, 20, 65, 90],
    },
    # theta=200: rho in [0, 1, 4, 16, 64]
    200: {
        "GENECONV": [0, 48, 88, 98, 100],
        "MAXCHI": [0, 45, 85, 98, 100],
        "CHIMAERA": [0, 48, 88, 98, 100],
        "HOMOPLASY": [0, 0, 0, 0, 12],
        "RETICULATE": [0, 35, 75, 95, 100],
        "RDP": [0, 32, 75, 96, 100],
        "SIMPLOT": [0, 5, 35, 60, 98],
    }
}

HISTORICAL_FP = {
    # alpha in [inf, 2.0, 0.5, 0.05]
    10: {
        "GENECONV": [0, 0, 0, 2],
        "MAXCHI": [0, 0, 0, 2],
        "CHIMAERA": [0, 0, 0, 2],
        "HOMOPLASY": [2, 0, 10, 28],
        "RETICULATE": [0, 0, 0, 2],
        "RDP": [0, 0, 0, 0],
    },
    50: {
        "GENECONV": [0, 0, 0, 4],
        "MAXCHI": [0, 0, 0, 4],
        "CHIMAERA": [0, 0, 0, 4],
        "HOMOPLASY": [0, 0, 5, 72],
        "RETICULATE": [0, 0, 2, 5],
        "RDP": [0, 0, 0, 2],
    },
    100: {
        "GENECONV": [0, 0, 2, 4],
        "MAXCHI": [0, 0, 2, 4],
        "CHIMAERA": [0, 0, 2, 4],
        "HOMOPLASY": [0, 2, 5, 86],
        "RETICULATE": [0, 2, 2, 5],
        "RDP": [0, 0, 0, 2],
    },
    200: {
        "GENECONV": [0, 0, 2, 5],
        "MAXCHI": [0, 0, 2, 5],
        "CHIMAERA": [0, 0, 2, 5],
        "HOMOPLASY": [0, 0, 2, 90],
        "RETICULATE": [0, 2, 5, 8],
        "RDP": [0, 0, 2, 5],
    }
}


def plot_posada_replication(summary_csv: Path, out_dir: Path):
    df = pd.read_csv(summary_csv)

    fig, axes = plt.subplots(4, 2, figsize=(11, 14), sharex="col", sharey=True)
    fig.subplots_adjust(hspace=0.22, wspace=0.15, top=0.93, bottom=0.12)

    thetas = [10.0, 50.0, 100.0, 200.0]
    rhos = [0.0, 1.0, 4.0, 16.0, 64.0]
    alphas_labels = [r"$\infty$", "2", "0.5", "0.05"]
    alphas_vals = ["inf", "2.0", "0.5", "0.05"]

    # Colors
    c_rhiz = "#D62728"   # Bold Crimson
    c_3seq = "#1F77B4"   # Royal Blue
    c_gene = "#2CA02C"   # Green
    c_max = "#FF7F0E"    # Orange
    c_homo = "#9467BD"   # Purple
    c_ret = "#8C564B"    # Brown
    c_rdp = "#17BECF"    # Teal

    for r_idx, theta in enumerate(thetas):
        # ---------------------------------------------------------------------
        # LEFT COLUMN: POWER
        # ---------------------------------------------------------------------
        ax_p = axes[r_idx, 0]
        ax_p.set_ylim(-2, 104)
        ax_p.grid(True, linestyle="--")

        # Slice summary data for power
        sub_p = df[(df["category"] == "power") & (df["theta"] == theta)].sort_values("rho")
        
        # Plot Historical Methods (faint lines)
        h_dict = HISTORICAL_POWER[int(theta)]
        ax_p.plot(range(5), h_dict["GENECONV"], color=c_gene, linestyle=":", marker="o", alpha=0.55, label="GENECONV")
        ax_p.plot(range(5), h_dict["MAXCHI"], color=c_max, linestyle=":", marker="s", alpha=0.55, label="MAXCHI")
        ax_p.plot(range(5), h_dict["HOMOPLASY"], color=c_homo, linestyle=":", marker="^", alpha=0.55, label="HOMOPLASY")
        ax_p.plot(range(5), h_dict["RETICULATE"], color=c_ret, linestyle=":", marker="v", alpha=0.55, label="RETICULATE")
        ax_p.plot(range(5), h_dict["RDP"], color=c_rdp, linestyle=":", marker="x", alpha=0.55, label="RDP")

        # Plot 3SEQ
        if "three_seq_detection_pct" in sub_p.columns:
            ax_p.plot(range(5), sub_p["three_seq_detection_pct"].values, color=c_3seq,
                      linestyle="-", marker="s", linewidth=2.0, label="3SEQ")

        # Plot RhizAeon (Primary)
        if "rhizaeon_detection_pct" in sub_p.columns:
            ax_p.plot(range(5), sub_p["rhizaeon_detection_pct"].values, color=c_rhiz,
                      linestyle="-", marker="D", linewidth=2.4, markersize=6, label="RhizAeon")

        # Panel Annotations
        div_approx = {10: "1\\%", 50: "5\\%", 100: "9\\%", 200: "17\\%"}[int(theta)]
        ax_p.text(0.04, 0.88, f"$\\theta = {int(theta)}$  (div $\\approx$ {div_approx})\n$\\alpha = \\infty$ (uniform rate)",
                  transform=ax_p.transAxes, fontsize=9.5, verticalalignment="top",
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#cccccc"))
        ax_p.set_ylabel("Detection Power (%)")

        # ---------------------------------------------------------------------
        # RIGHT COLUMN: FALSE POSITIVES (rho = 0)
        # ---------------------------------------------------------------------
        ax_fp = axes[r_idx, 1]
        ax_fp.set_ylim(-2, 104)
        ax_fp.grid(True, linestyle="--")

        # Historical FP
        h_fp = HISTORICAL_FP[int(theta)]
        ax_fp.plot(range(4), h_fp["HOMOPLASY"], color=c_homo, linestyle=":", marker="^", alpha=0.6, label="HOMOPLASY")
        ax_fp.plot(range(4), h_fp["GENECONV"], color=c_gene, linestyle=":", marker="o", alpha=0.55, label="GENECONV")
        ax_fp.plot(range(4), h_fp["MAXCHI"], color=c_max, linestyle=":", marker="s", alpha=0.55, label="MAXCHI")
        ax_fp.plot(range(4), h_fp["RETICULATE"], color=c_ret, linestyle=":", marker="v", alpha=0.55, label="RETICULATE")
        ax_fp.plot(range(4), h_fp["RDP"], color=c_rdp, linestyle=":", marker="x", alpha=0.55, label="RDP")

        # Get summary data for FP: alpha in [inf, 2.0, 0.5, 0.05]
        fp_rhiz = []
        fp_3seq = []
        for a_val in alphas_vals:
            if a_val == "inf":
                row = df[(df["category"] == "power") & (df["theta"] == theta) & (df["rho"] == 0.0)]
            else:
                row = df[(df["category"] == "false_positive") & (df["theta"] == theta) & (df["alpha"] == float(a_val))]
            if len(row) > 0:
                fp_rhiz.append(row["rhizaeon_detection_pct"].values[0])
                fp_3seq.append(row["three_seq_detection_pct"].values[0])
            else:
                fp_rhiz.append(0.0)
                fp_3seq.append(0.0)

        ax_fp.plot(range(4), fp_3seq, color=c_3seq, linestyle="-", marker="s", linewidth=2.0, label="3SEQ")
        ax_fp.plot(range(4), fp_rhiz, color=c_rhiz, linestyle="-", marker="D", linewidth=2.4, markersize=6, label="RhizAeon")

        # 5% nominal threshold reference line
        ax_fp.axhline(5.0, color="#888888", linestyle="--", linewidth=1.0, alpha=0.7)

        ax_fp.text(0.04, 0.88, f"$\\theta = {int(theta)}$\n$\\rho = 0$ (no recombination)",
                  transform=ax_fp.transAxes, fontsize=9.5, verticalalignment="top",
                  bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="#cccccc"))
        ax_fp.set_ylabel("False Positive Rate (%)")

    # Set X tick labels for bottom row
    axes[3, 0].set_xticks(range(5))
    axes[3, 0].set_xticklabels(["0", "1", "4", "16", "64"])
    axes[3, 0].set_xlabel("Recombination parameter ($\\rho = 4Nrl$)")

    axes[3, 1].set_xticks(range(4))
    axes[3, 1].set_xticklabels(alphas_labels)
    axes[3, 1].set_xlabel("Rate variation among sites ($\\alpha$)")

    # Column titles
    axes[0, 0].set_title("POWER (Simulations I)", fontsize=12, fontweight="bold", pad=12)
    axes[0, 1].set_title("FALSE POSITIVES (Simulations II)", fontsize=12, fontweight="bold", pad=12)

    # Master figure title
    fig.suptitle("Posada & Crandall (2001) PNAS Benchmark Replication:\nRhizAeon and 3SEQ vs Historical Recombination Detection Methods",
                 fontsize=13, fontweight="bold", y=0.98)

    # Master legend below plots
    handles, labels = axes[0, 0].get_legend_handles_labels()
    # Deduplicate
    by_label = dict(zip(labels, handles))
    # Order: RhizAeon, 3SEQ, then historical
    order = ["RhizAeon", "3SEQ", "GENECONV", "MAXCHI", "HOMOPLASY", "RETICULATE", "RDP"]
    ordered_handles = [by_label[k] for k in order if k in by_label]
    ordered_labels = [k for k in order if k in by_label]

    fig.legend(ordered_handles, ordered_labels, loc="lower center", ncol=7,
               bbox_to_anchor=(0.5, 0.03), frameon=True, facecolor="#f9f9f9", edgecolor="#cccccc", fontsize=9.5)

    png_path = out_dir / "fig_posada_2001_replication.png"
    pdf_path = out_dir / "fig_posada_2001_replication.pdf"
    plt.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[*] Generated replication figure: {png_path} and {pdf_path}")


if __name__ == "__main__":
    out_dir = Path("benchmarks/posada_2001_benchmark")
    summary_csv = out_dir / "posada_summary.csv"
    if summary_csv.exists():
        plot_posada_replication(summary_csv, out_dir)
    else:
        print(f"Summary CSV not found at: {summary_csv}")
