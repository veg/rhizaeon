"""
simulations/02_rivet_sars_cov_2/plot_rivet_benchmark.py
======================================================
BioVis-compliant multi-panel figure generator for the RIVET SARS-CoV-2 benchmark.
Generates paper/Figures/fig_rivet_benchmark.pdf (and PNG).
"""

import os
import re
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STUDY_DIR = Path(__file__).resolve().parent
FIGURES_DIR = PROJECT_ROOT / "paper" / "Figures"
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# BioVis Palette
COLOR_RHIZ = "#1f77b4"     # RhizAeon Blue
COLOR_3SEQ = "#ff7f0e"     # 3SEQ Orange
COLOR_RIVET = "#2ca02c"    # RIVET Green
COLOR_GRAY = "#7f7f7f"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.titlesize": 11,
    "axes.linewidth": 0.8,
    "grid.linewidth": 0.5,
    "grid.alpha": 0.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


def plot_rivet_figure():
    summary_path = STUDY_DIR / "rivet_summary.csv"
    emp_path = STUDY_DIR / "rivet_empirical_trios.csv"
    raw_path = STUDY_DIR / "rivet_raw_results.csv"

    if not summary_path.exists() or not emp_path.exists():
        print(f"Error: Missing data files in {STUDY_DIR}")
        return

    df_sum = pd.read_csv(summary_path)
    df_emp = pd.read_csv(emp_path)
    df_raw = pd.read_csv(raw_path)

    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.8))
    plt.subplots_adjust(wspace=0.30, hspace=0.34, left=0.11, right=0.96, top=0.93, bottom=0.08)

    # -------------------------------------------------------------
    # Panel A: Statistical Power vs Parental Divergence (1-BP)
    # -------------------------------------------------------------
    ax_a = axes[0, 0]
    df_1bp = df_sum[df_sum["n_bps"] == 1]
    d_vals = sorted(df_1bp["d"].unique())

    # Aggregate over m for main curves, and plot m=0 vs m=3 as shaded envelope
    mean_rhiz = []
    mean_3seq = []
    rhiz_m0 = []
    rhiz_m3 = []
    three_m0 = []
    three_m3 = []

    for d in d_vals:
        sub = df_1bp[df_1bp["d"] == d]
        mean_rhiz.append(sub["rhizaeon_power"].mean() * 100)
        mean_3seq.append(sub["three_seq_power"].mean() * 100)
        rhiz_m0.append(sub[sub["m"] == 0]["rhizaeon_power"].values[0] * 100)
        rhiz_m3.append(sub[sub["m"] == 3]["rhizaeon_power"].values[0] * 100)
        three_m0.append(sub[sub["m"] == 0]["three_seq_power"].values[0] * 100)
        three_m3.append(sub[sub["m"] == 3]["three_seq_power"].values[0] * 100)

    ax_a.plot(d_vals, mean_rhiz, "o-", color=COLOR_RHIZ, lw=1.8, ms=5, label="RhizAeon (Mean)", zorder=4)
    ax_a.fill_between(d_vals, rhiz_m3, rhiz_m0, color=COLOR_RHIZ, alpha=0.15, label="RhizAeon ($m=0..3$)")

    ax_a.plot(d_vals, mean_3seq, "s--", color=COLOR_3SEQ, lw=1.8, ms=5, label="3SEQ (Mean)", zorder=3)
    ax_a.fill_between(d_vals, three_m3, three_m0, color=COLOR_3SEQ, alpha=0.15, label="3SEQ ($m=0..3$)")

    ax_a.set_title("A  Recombination Power vs Parental Divergence", fontweight="semibold", loc="left")
    ax_a.set_xlabel("Parental Mutational Divergence $d$ (SNPs)")
    ax_a.set_ylabel("Detection Power (%)")
    ax_a.set_ylim(-2, 105)
    ax_a.set_xticks(d_vals)
    ax_a.grid(True, linestyle=":")
    ax_a.legend(loc="lower right", frameon=True, framealpha=0.9, edgecolor="none", fontsize=7)

    # -------------------------------------------------------------
    # Panel B: Breakpoint Localization Error (Spatial Precision)
    # -------------------------------------------------------------
    ax_b = axes[0, 1]
    err_data = []
    pos_labels = []
    for d in d_vals:
        sub_raw = df_raw[(df_raw["n_bps"] == 1) & (df_raw["d"] == d) & (df_raw["rhizaeon_detected"])]
        errs = sub_raw["rhizaeon_spatial_err_nt"].dropna().values
        if len(errs) > 0:
            err_data.append(errs)
            pos_labels.append(f"$d={d}$")

    bp = ax_b.boxplot(err_data, positions=range(len(pos_labels)), widths=0.45,
                      patch_artist=True, showfliers=False,
                      medianprops=dict(color="black", lw=1.5),
                      boxprops=dict(facecolor=COLOR_RHIZ, alpha=0.6, color=COLOR_RHIZ))

    ax_b.set_title("B  Spatial Breakpoint Error Across 29.9 kb", fontweight="semibold", loc="left")
    ax_b.set_xlabel("Divergence Scenario")
    ax_b.set_ylabel("Spatial Error $|\hat{b} - b_{\mathrm{true}}|$ (nt)")
    ax_b.set_xticks(range(len(pos_labels)))
    ax_b.set_xticklabels(pos_labels)
    ax_b.grid(True, linestyle=":", axis="y")

    # Annotate mean errors
    for idx, d in enumerate(d_vals):
        sub_raw = df_raw[(df_raw["n_bps"] == 1) & (df_raw["d"] == d) & (df_raw["rhizaeon_detected"])]
        m_err = sub_raw["rhizaeon_spatial_err_nt"].dropna().mean()
        if not np.isnan(m_err):
            ax_b.text(idx, m_err + 80, f"{m_err:.0f} nt", ha="center", fontsize=7, color="#0c3c60", fontweight="bold")

    # -------------------------------------------------------------
    # Panel C: Empirical Pandemic Concordance (481 PASS Trios)
    # -------------------------------------------------------------
    ax_c = axes[1, 0]
    # Extract RIVET breakpoint midpoints and RhizAeon detected midpoints
    rivet_bps = []
    rhiz_emp_bps = []
    for _, r in df_emp.iterrows():
        m = re.search(r"\((\d+)\s*,\s*(\d+)\)", str(r["rivet_bp1_interval"]))
        if m:
            l, rt = int(m.group(1)), int(m.group(2))
            if l > 0:
                rivet_bps.append((l + rt) / 2.0)
        if r["rhizaeon_detected"]:
            for c in str(r["rhizaeon_bps"]).split(";"):
                if c.strip():
                    try:
                        rhiz_emp_bps.append(float(c))
                    except ValueError:
                        pass

    bins = np.linspace(0, 29903, 35)
    ax_c.hist(rivet_bps, bins=bins, color=COLOR_RIVET, alpha=0.45, label=f"RIVET Intervals ($N={len(rivet_bps)}$)", density=True)
    ax_c.hist(rhiz_emp_bps, bins=bins, color=COLOR_RHIZ, alpha=0.55, label=f"RhizAeon Detected ($N={len(rhiz_emp_bps)}$)", density=True, histtype="step", lw=1.8)

    # Spike gene highlight (21563 - 25384)
    ax_c.axvspan(21563, 25384, color="#d95f02", alpha=0.15, label="Spike Gene ($S$)")

    ax_c.set_title("C  Empirical Genomic Breakpoint Distribution", fontweight="semibold", loc="left")
    ax_c.set_xlabel("Genomic Position (nt, Wuhan-Hu-1 Coordinates)")
    ax_c.set_ylabel("Breakpoint Density")
    ax_c.set_xlim(0, 29903)
    ax_c.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, p: f"{int(x/1000)}k"))
    ax_c.grid(True, linestyle=":")
    ax_c.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="none", fontsize=7)

    # -------------------------------------------------------------
    # Panel D: Computational Latency Scaling (Whole-Genome Speed)
    # -------------------------------------------------------------
    ax_d = axes[1, 1]
    methods = ["RhizAeon", "3SEQ"]
    emp_lat = [df_emp["rhizaeon_runtime_ms"].mean(), df_emp["three_seq_runtime_ms"].mean()]
    sim_lat = [df_raw["rhizaeon_runtime_ms"].mean(), df_raw["three_seq_runtime_ms"].mean()]

    x = np.arange(len(methods))
    width = 0.35

    rects1 = ax_d.bar(x - width/2, emp_lat, width, label="Empirical Trios ($N=481$)", color=["#1f77b4", "#ff7f0e"], alpha=0.85)
    rects2 = ax_d.bar(x + width/2, sim_lat, width, label="Simulated Trios ($N=1900$)", color=["#aec7e8", "#ffbb78"], edgecolor=["#1f77b4", "#ff7f0e"], lw=1.2)

    ax_d.set_title("D  Computational Latency per 30 kb Genome", fontweight="semibold", loc="left")
    ax_d.set_ylabel("Execution Latency (ms / genome)")
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(methods, fontweight="semibold")
    ax_d.set_yscale("log")
    ax_d.set_ylim(1, 150)
    ax_d.grid(True, linestyle=":", axis="y")
    ax_d.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="none", fontsize=7)

    # Annotate speedup
    speedup = sim_lat[1] / sim_lat[0]
    ax_d.text(0.5, 0.72, f"RhizAeon Speedup:\n{speedup:.1f}x Faster than 3SEQ", transform=ax_d.transAxes,
              ha="center", va="center", bbox=dict(boxstyle="round,pad=0.3", fc="#eef4f8", ec=COLOR_RHIZ, lw=0.8),
              fontsize=7.5, fontweight="semibold", color="#0c3c60")

    # Save figure to paper/Figures and study directory
    out_pdf = FIGURES_DIR / "fig_rivet_benchmark.pdf"
    out_png = FIGURES_DIR / "fig_rivet_benchmark.png"
    plt.savefig(out_pdf, dpi=300)
    plt.savefig(out_png, dpi=300)

    # Copy to study directory
    plt.savefig(STUDY_DIR / "fig_rivet_benchmark.pdf", dpi=300)
    plt.savefig(STUDY_DIR / "fig_rivet_benchmark.png", dpi=300)
    plt.close()

    print(f"Generated BioVis-compliant figure: {out_pdf} and {out_png}")


if __name__ == "__main__":
    plot_rivet_figure()
