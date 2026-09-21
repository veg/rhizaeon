#!/usr/bin/env python3
"""
benchmarks/recombinhunt_2024/plot_recombinhunt_benchmark.py
===========================================================
Generates publication-quality figure comparing RhizAeon vs RecombinHunt
across the 10,500-genome SARS-CoV-2 viral benchmark suite (Alfonsi et al. 2024).
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
SUMMARY_CSV = SCRIPT_DIR / "rhizaeon_recombinhunt_summary.csv"
RAW_CSV = SCRIPT_DIR / "rhizaeon_recombinhunt_raw_results.csv"
FIG_PNG = REPO_ROOT / "paper" / "Figures" / "fig_recombinhunt_benchmark.png"
FIG_PDF = REPO_ROOT / "paper" / "Figures" / "fig_recombinhunt_benchmark.pdf"

# Set publication style
plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10.5,
    "ytick.labelsize": 10.5,
    "legend.fontsize": 10.5,
    "figure.titlesize": 15,
    "lines.linewidth": 2.2,
    "lines.markersize": 7.0,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})

def plot_benchmark():
    sum_df = pd.read_csv(SUMMARY_CSV)
    raw_df = pd.read_csv(RAW_CSV)

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.28)

    color_rhiz = "#1b9e77"    # Teal / Green
    color_rh = "#d95f02"      # Dark Orange
    color_3seq = "#7570b3"    # Purple

    # -------------------------------------------------------------
    # PANEL A: 1-BP Sensitivity vs Noise Level
    # -------------------------------------------------------------
    ax1 = fig.add_subplot(gs[0, 0])
    df_1bp = sum_df[sum_df['scenario_type'] == '1BP_Recombinant'].sort_values('noise_level')
    
    ax1.plot(df_1bp['noise_level'], df_1bp['rhizaeon_sensitivity_pct'], marker='o', color=color_rhiz, label='RhizAeon (Prefix FDA)')
    ax1.plot(df_1bp['noise_level'], df_1bp['recombinhunt_sensitivity_pct'], marker='s', linestyle='--', color=color_rh, label='RecombinHunt (Alfonsi 2024)')
    
    ax1.set_title("(A) 1-BP Sensitivity vs. Noise Level", fontweight="bold")
    ax1.set_xlabel("Added Noise Mutations (Count)")
    ax1.set_ylabel("Detection Sensitivity (%)")
    ax1.set_ylim(85, 102)
    ax1.set_xticks(df_1bp['noise_level'])
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend(loc="lower left", frameon=True)

    # -------------------------------------------------------------
    # PANEL B: 2-BP Sensitivity vs Noise Level
    # -------------------------------------------------------------
    ax2 = fig.add_subplot(gs[0, 1])
    df_2bp = sum_df[sum_df['scenario_type'] == '2BP_Recombinant'].sort_values('noise_level')

    ax2.plot(df_2bp['noise_level'], df_2bp['rhizaeon_sensitivity_pct'], marker='o', color=color_rhiz, label='RhizAeon (Prefix FDA)')
    ax2.plot(df_2bp['noise_level'], df_2bp['recombinhunt_sensitivity_pct'], marker='s', linestyle='--', color=color_rh, label='RecombinHunt (Alfonsi 2024)')

    ax2.set_title("(B) 2-BP Sensitivity vs. Noise Level", fontweight="bold")
    ax2.set_xlabel("Added Noise Mutations (Count)")
    ax2.set_ylabel("Detection Sensitivity (%)")
    ax2.set_ylim(70, 102)
    ax2.set_xticks(df_2bp['noise_level'])
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend(loc="lower left", frameon=True)

    # -------------------------------------------------------------
    # PANEL C: False Positive Rate vs Noise Level (0-BP Null)
    # -------------------------------------------------------------
    ax3 = fig.add_subplot(gs[0, 2])
    df_0bp = sum_df[sum_df['scenario_type'] == '0BP_Null'].sort_values('noise_level')

    ax3.plot(df_0bp['noise_level'], df_0bp['rhizaeon_fpr_pct'], marker='o', color=color_rhiz, label='RhizAeon (FPR: 0.0%)')
    ax3.plot(df_0bp['noise_level'], df_0bp['recombinhunt_fpr_pct'], marker='s', linestyle='--', color=color_rh, label='RecombinHunt')

    ax3.set_title("(C) False Positive Rate vs. Noise (0-BP Null)", fontweight="bold")
    ax3.set_xlabel("Added Noise Mutations (Count)")
    ax3.set_ylabel("False Positive Rate (%)")
    ax3.set_ylim(-0.5, 10.0)
    ax3.set_xticks(df_0bp['noise_level'])
    ax3.grid(True, linestyle=':', alpha=0.6)
    ax3.legend(loc="upper left", frameon=True)

    # -------------------------------------------------------------
    # PANEL D: Breakpoint Spatial Precision (1-BP & 2-BP)
    # -------------------------------------------------------------
    ax4 = fig.add_subplot(gs[1, 0])
    noise_levels = sorted(df_1bp['noise_level'].unique())
    spatial_1bp = df_1bp['rhizaeon_spatial_accuracy_pct'].values
    spatial_2bp = df_2bp['rhizaeon_spatial_accuracy_pct'].values

    x = np.arange(len(noise_levels))
    width = 0.35

    ax4.bar(x - width/2, spatial_1bp, width, label='RhizAeon 1-BP', color=color_rhiz, alpha=0.85)
    ax4.bar(x + width/2, spatial_2bp, width, label='RhizAeon 2-BP (Both BPs)', color="#2b83ba", alpha=0.85)

    ax4.axhline(99.4, color=color_rh, linestyle='--', label='RecombinHunt 1-BP Avg (99.4%)')
    ax4.axhline(97.6, color="#e7298a", linestyle=':', label='RecombinHunt 2-BP Avg (97.6%)')

    ax4.set_title("(D) Spatial Accuracy (Exact Ground Truth Interval)", fontweight="bold")
    ax4.set_xlabel("Added Noise Mutations (Count)")
    ax4.set_ylabel("In Ground Truth Interval (%)")
    ax4.set_ylim(85, 105)
    ax4.set_xticks(x)
    ax4.set_xticklabels(noise_levels)
    ax4.grid(True, linestyle=':', alpha=0.6)
    ax4.legend(loc="lower left", frameon=True, fontsize=9.5)

    # -------------------------------------------------------------
    # PANEL E: Spatial Error Distribution (2-BP Challenging Cases)
    # -------------------------------------------------------------
    ax5 = fig.add_subplot(gs[1, 1])
    raw_2bp = raw_df[raw_df['true_bp_num'] == 2]
    errs_bp1 = raw_2bp['spatial_err_bp1'].dropna().values
    errs_bp2 = raw_2bp['spatial_err_bp2'].dropna().values
    all_errs = np.concatenate([errs_bp1, errs_bp2])

    ax5.hist(all_errs[all_errs > 0], bins=30, color="#2b83ba", edgecolor="black", alpha=0.75)
    ax5.axvline(0, color=color_rhiz, linewidth=3, label=f'Exact Interval: {(all_errs == 0).mean()*100:.1f}%')
    ax5.set_title("(E) Spatial Error Distribution on Divergent Cases", fontweight="bold")
    ax5.set_xlabel("Spatial Distance from True Interval (nt)")
    ax5.set_ylabel("Breakpoint Count")
    ax5.grid(True, linestyle=':', alpha=0.6)
    ax5.legend(loc="upper right", frameon=True)

    # -------------------------------------------------------------
    # PANEL F: Throughput and Scalability Comparison
    # -------------------------------------------------------------
    ax6 = fig.add_subplot(gs[1, 2])
    methods = ['RhizAeon', 'RecombinHunt', '3SEQ (Est.)', 'GARD (Est.)']
    # RhizAeon: 1.27 ms per genome
    # RecombinHunt: ~45 ms per genome (approx 8 minutes for 10.5k reported)
    # 3SEQ: ~25 ms per triplet
    # GARD: ~300,000 ms (5 minutes) per genome
    times_ms = [1.27, 45.0, 25.0, 300000.0]
    colors_f = [color_rhiz, color_rh, color_3seq, "#e41a1c"]

    bars = ax6.bar(methods, times_ms, color=colors_f, alpha=0.85, edgecolor="black")
    ax6.set_yscale('log')
    ax6.set_title("(F) Computational Scalability (Runtime / Genome)", fontweight="bold")
    ax6.set_ylabel("Runtime per 30-kb Genome (ms, log scale)")
    ax6.grid(True, linestyle=':', alpha=0.6)

    # Annotate bars
    for bar, val in zip(bars, times_ms):
        h = bar.get_height()
        ax6.annotate(f'{val:.2f} ms' if val < 1000 else f'{val/1000:.0f} s',
                     xy=(bar.get_x() + bar.get_width() / 2, h),
                     xytext=(0, 4), textcoords="offset points",
                     ha='center', va='bottom', fontsize=9.5, fontweight='bold')

    fig.suptitle("RhizAeon vs. RecombinHunt: 10,500-Genome SARS-CoV-2 Benchmark Replication\n(Alfonsi et al. 2024, Nature Communications 15:3313)",
                 fontsize=15, fontweight="bold", y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(FIG_PDF)
    plt.close()
    
    # Render PNG using sips from PDF for crisp rendering without Agg driver mismatch
    import subprocess
    subprocess.run(["/usr/bin/sips", "-s", "format", "png", str(FIG_PDF), "--out", str(FIG_PNG)], check=True)
    
    # Also save copies in local simulation directory
    local_pdf = SCRIPT_DIR / "fig_recombinhunt_benchmark.pdf"
    local_png = SCRIPT_DIR / "fig_recombinhunt_benchmark.png"
    subprocess.run(["cp", str(FIG_PDF), str(local_pdf)], check=True)
    subprocess.run(["cp", str(FIG_PNG), str(local_png)], check=True)
    print(f"Saved publication figures to:\n  {FIG_PNG}\n  {FIG_PDF}\n  {local_png}\n  {local_pdf}")

if __name__ == "__main__":
    plot_benchmark()
