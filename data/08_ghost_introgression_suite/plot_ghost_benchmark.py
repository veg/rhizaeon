#!/usr/bin/env python3
"""
simulations/08_ghost_introgression_suite/plot_ghost_benchmark.py
================================================================
BioVis publication-quality multipanel figure evaluating RhizAeon performance
under incomplete parental lineage sampling, donor drift, and deep taxonomic scaling:
- Panel A: Event Detection Power vs Donor Drift (Delta d in [0%, 10%]) and Complete Ghost.
- Panel B: Dual-Breakpoint Recovery and Spatial Error (MAE in nt) vs Donor Drift.
- Panel C: Standardized Manifold Procrustes Kinetic Z-Score Surge vs Donor Drift.
- Panel D: Computational Wall-Clock Scaling vs Cohort Size N (5 to 100 taxa).
"""

import sys
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# BioVis font and styling configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
    "lines.linewidth": 2.2,
    "lines.markersize": 7.0,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
GHOST_CSV = SCRIPT_DIR / "ghost_summary.csv"
SCALING_CSV = SCRIPT_DIR / "scaling_summary.csv"
FIG_PDF = REPO_ROOT / "paper" / "Figures" / "fig_ghost_benchmark.pdf"
FIG_PNG = REPO_ROOT / "paper" / "Figures" / "fig_ghost_benchmark.png"

# Okabe-Ito Color Palette
C_PRIMARY = "#0284c7"   # Blue
C_VERMILION = "#d95f02" # Vermilion
C_EMERALD = "#059669"   # Emerald green
C_PURPLE = "#7570b3"    # Purple
C_GRAY = "#64748b"      # Gray

def plot_ghost_benchmark():
    ghost_df = pd.read_csv(GHOST_CSV)
    scaling_df = pd.read_csv(SCALING_CSV)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))
    plt.subplots_adjust(hspace=0.32, wspace=0.28)

    # Filter non-complete ghost for continuous drift curve
    drift_df = ghost_df[ghost_df["is_complete_ghost"] == False].sort_values("donor_divergence")
    drift_pct = drift_df["donor_divergence"] * 100.0

    # -----------------------------------------------------------------
    # Panel A: Event Detection Power vs Donor Drift
    # -----------------------------------------------------------------
    ax1 = axes[0, 0]
    ax1.plot(drift_pct, drift_df["event_power_pct"], marker="o", color=C_PRIMARY, label="Event Detection Power")
    # Mark complete ghost point
    ghost_row = ghost_df[ghost_df["is_complete_ghost"] == True]
    if not ghost_row.empty:
        ax1.scatter([12.0], [ghost_row["event_power_pct"].values[0]], color=C_VERMILION, s=90, zorder=5, marker="X",
                    label=f"100% Missing Clade: {ghost_row['event_power_pct'].values[0]:.0f}%")

    ax1.set_title("A. Detection Power Under Donor Drift", loc="left")
    ax1.set_xlabel("Parental Donor Drift $\\Delta d$ (% Divergence)")
    ax1.set_ylabel("Detection Power (%)")
    ax1.set_ylim(-5, 110)
    ax1.set_xlim(-0.5, 13.5)
    ax1.set_xticks([0, 1, 3, 5, 10, 12])
    ax1.set_xticklabels(["0%", "1%", "3%", "5%", "10%", "Excised"])
    ax1.grid(True, linestyle="--")
    ax1.legend(frameon=True, loc="lower left")

    # -----------------------------------------------------------------
    # Panel B: Dual-Breakpoint Power & Spatial Accuracy
    # -----------------------------------------------------------------
    ax2 = axes[0, 1]
    ax2.plot(drift_pct, drift_df["dual_bp_power_pct"], marker="s", color=C_EMERALD, label="Dual-BP Recovery (%)")
    ax2.set_xlabel("Parental Donor Drift $\\Delta d$ (% Divergence)")
    ax2.set_ylabel("Dual-Breakpoint Recovery (%)", color=C_EMERALD)
    ax2.tick_params(axis="y", labelcolor=C_EMERALD)
    ax2.set_ylim(-5, 105)
    ax2.set_xlim(-0.5, 11)
    ax2.grid(True, linestyle="--")

    ax2_twin = ax2.twinx()
    ax2_twin.plot(drift_pct, drift_df["breakpoint_mae_nt"], marker="^", linestyle="--", color=C_VERMILION, label="Spatial Error (MAE nt)")
    ax2_twin.set_ylabel("Spatial Error (MAE in nt)", color=C_VERMILION)
    ax2_twin.tick_params(axis="y", labelcolor=C_VERMILION)
    ax2_twin.set_ylim(0, 30)

    # Combined legend
    lines1, labels1 = ax2.get_legend_handles_labels()
    lines2, labels2 = ax2_twin.get_legend_handles_labels()
    ax2.legend(lines1 + lines2, labels1 + labels2, frameon=True, loc="center right")
    ax2.set_title("B. Dual-Breakpoint Resolution and Spatial Precision", loc="left")

    # -----------------------------------------------------------------
    # Panel C: Manifold Procrustes Kinetic Z-Score Surge
    # -----------------------------------------------------------------
    ax3 = axes[1, 0]
    ax3.plot(drift_pct, drift_df["max_z"], marker="D", color=C_PURPLE, lw=2.4, markersize=8, label="Max Kinetic $Z$-Score")
    ax3.axhline(3.0, color=C_GRAY, linestyle=":", lw=1.5, label="Statistical Threshold ($Z = 3.0$)")
    ax3.set_title("C. Manifold Procrustes Kinetic Residual Surge", loc="left")
    ax3.set_xlabel("Parental Donor Drift $\\Delta d$ (% Divergence)")
    ax3.set_ylabel("Standardized Kinetic $Z$-Score")
    ax3.set_ylim(0, 10.0)
    ax3.set_xlim(-0.5, 11)
    ax3.grid(True, linestyle="--")
    ax3.legend(frameon=True, loc="upper left")

    # -----------------------------------------------------------------
    # Panel D: Computational Cohort Scaling vs N
    # -----------------------------------------------------------------
    ax4 = axes[1, 1]
    scale_5pct = scaling_df[scaling_df["divergence"] == 0.05].sort_values("N")
    scale_2pct = scaling_df[scaling_df["divergence"] == 0.02].sort_values("N")

    ax4.plot(scale_5pct["N"], scale_5pct["runtime_ms"], marker="o", color=C_PRIMARY, label="RhizAeon ($d = 0.05$)")
    ax4.plot(scale_2pct["N"], scale_2pct["runtime_ms"], marker="^", linestyle="--", color=C_GRAY, label="RhizAeon ($d = 0.02$)")

    # Triplet combinatorial scaling reference: ~ O(N^3)
    N_ref = np.array([5, 10, 20, 50, 100])
    triplet_ref = (N_ref ** 3) / 1000.0 * 1.5  # ms
    ax4.plot(N_ref, triplet_ref, ":", color=C_VERMILION, lw=1.8, label="Triplet Enumeration $\\mathcal{O}(N^3)$ (Est.)")

    ax4.set_yscale("log")
    ax4.set_title("D. Computational Scaling Across Deep Cohorts", loc="left")
    ax4.set_xlabel("Taxonomic Cohort Size ($N$ Isolates)")
    ax4.set_ylabel("Execution Time per Alignment (ms, log scale)")
    ax4.grid(True, linestyle="--", which="both")
    ax4.legend(frameon=True, loc="upper left")

    # Save outputs
    plt.savefig(FIG_PDF, bbox_inches="tight")
    plt.close()

    # Render high-resolution PNG via sips
    subprocess.run(["/usr/bin/sips", "-s", "format", "png", str(FIG_PDF), "--out", str(FIG_PNG)], check=True)

    # Save local copies
    local_pdf = SCRIPT_DIR / "fig_ghost_benchmark.pdf"
    local_png = SCRIPT_DIR / "fig_ghost_benchmark.png"
    subprocess.run(["cp", str(FIG_PDF), str(local_pdf)], check=True)
    subprocess.run(["cp", str(FIG_PNG), str(local_png)], check=True)
    print(f"Generated Figures:\n  {FIG_PDF}\n  {FIG_PNG}\n  {local_pdf}\n  {local_png}")

if __name__ == "__main__":
    plot_ghost_benchmark()
