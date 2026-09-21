"""
simulations/02_phipack_bruen_2006/plot_phipack_benchmark.py
===========================================================
BioVis publication-quality multipanel figure replicating Bruen et al. (2006):
- Panel A: Detection Power vs Recombination Intensity rho (n=10 vs n=50).
- Panel B: False Positive Rates across Null Regimes (Constant, Growth, Gamma).
- Panel C: Head-to-head Power across Methods (RhizAeon, PHI, MaxChi, NSS, 3SEQ).
- Panel D: Execution Latency (ms) comparison.
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

# Color palette (BioVis colorblind-friendly)
COLORS = {
    "RhizAeon": "#0284c7",   # Sky blue
    "PHI": "#059669",        # Emerald green
    "MaxChi": "#d97706",     # Amber
    "NSS": "#7c3aed",        # Purple
    "3SEQ": "#dc2626"        # Rose red
}


def plot_phipack_figure(summary_csv: Path, raw_csv: Path, out_pdf: Path, out_png: Path):
    df_sum = pd.read_csv(summary_csv)
    df_raw = pd.read_csv(raw_csv)

    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.8), dpi=300)
    plt.subplots_adjust(wspace=0.28, hspace=0.35, top=0.92, bottom=0.08, left=0.08, right=0.96)

    # -------------------------------------------------------------
    # Panel A: Detection Power vs rho (n=10 vs n=50, theta=20)
    # -------------------------------------------------------------
    ax_a = axes[0, 0]
    p_df = df_sum[(df_sum["category"] == "power") & (df_sum["theta"] == 20.0)]
    
    # n = 50 curves
    p50 = p_df[p_df["n_taxa"] == 50].sort_values("rho")
    p10 = p_df[p_df["n_taxa"] == 10].sort_values("rho")

    method_map = {
        "RhizAeon": "rhiz_rate",
        "PHI": "phi_rate",
        "MaxChi": "maxchi_rate",
        "NSS": "nss_rate",
        "3SEQ": "three_rate"
    }

    for method, col, marker in [
        ("RhizAeon", COLORS["RhizAeon"], "o"),
        ("PHI", COLORS["PHI"], "s"),
        ("MaxChi", COLORS["MaxChi"], "^"),
        ("NSS", COLORS["NSS"], "d"),
        ("3SEQ", COLORS["3SEQ"], "v")
    ]:
        rate_col = method_map[method]
        ax_a.plot(p50["rho"], p50[rate_col], label=f"{method} (n=50)", color=col, marker=marker, markersize=4.5)
        ax_a.plot(p10["rho"], p10[rate_col], linestyle="--", alpha=0.55, color=col, marker=marker, markersize=3.5)

    ax_a.set_xscale("log", base=2)
    ax_a.set_xticks([1, 4, 16, 64])
    ax_a.get_xaxis().set_major_formatter(ticker.ScalarFormatter())
    ax_a.set_xlabel(r"Population Recombination Rate ($\rho = 4N_0 r$)")
    ax_a.set_ylabel("Detection Power (%)")
    ax_a.set_title("A   Detection Power vs Recombination Rate", loc="left", fontweight="bold")
    ax_a.set_ylim(-2, 105)
    ax_a.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1")
    ax_a.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0", ncol=2, fontsize=7.0)

    # -------------------------------------------------------------
    # Panel B: False Positive Rates across Null Regimes
    # -------------------------------------------------------------
    ax_b = axes[0, 1]
    null_cats = [
        ("Neutral Constant", df_raw[df_raw["category"] == "null_constant"]),
        (r"Exp. Growth ($\beta=5$)", df_raw[(df_raw["category"] == "null_growth") & (df_raw["beta"] == 5.0)]),
        (r"Exp. Growth ($\beta=20$)", df_raw[(df_raw["category"] == "null_growth") & (df_raw["beta"] == 20.0)]),
        (r"Gamma ($\alpha=2.0$)", df_raw[(df_raw["category"] == "null_gamma") & (df_raw["alpha"] == 2.0)]),
        (r"Gamma ($\alpha=0.5$)", df_raw[(df_raw["category"] == "null_gamma") & (df_raw["alpha"] == 0.5)])
    ]

    labels = [c[0] for c in null_cats]
    x_indices = np.arange(len(labels))
    width = 0.16

    methods = [("RhizAeon", "rhiz_detected", COLORS["RhizAeon"]),
               ("PHI", "phi_detected", COLORS["PHI"]),
               ("MaxChi", "maxchi_detected", COLORS["MaxChi"]),
               ("NSS", "nss_detected", COLORS["NSS"]),
               ("3SEQ", "three_detected", COLORS["3SEQ"])]

    for idx, (m_name, m_col, col) in enumerate(methods):
        fpr_vals = [c[1][m_col].mean() * 100.0 for c in null_cats]
        ax_b.bar(x_indices + (idx - 2) * width, fpr_vals, width, label=m_name, color=col, alpha=0.88, edgecolor="none")

    ax_b.axhline(5.0, color="#dc2626", linestyle="--", linewidth=1.0, alpha=0.7, label=r"Nominal $\alpha = 0.05$")
    ax_b.set_xticks(x_indices)
    ax_b.set_xticklabels(labels, rotation=25, ha="right", fontsize=7.5)
    ax_b.set_ylabel("False Positive Rate (%)")
    ax_b.set_title("B   False Positive Calibration under Null Regimes", loc="left", fontweight="bold")
    ax_b.set_ylim(0, 30)
    ax_b.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_b.legend(loc="upper right", frameon=True, framealpha=0.9, edgecolor="#e2e8f0", fontsize=7.0, ncol=2)

    # -------------------------------------------------------------
    # Panel C: Method Sensitivity at High Recombination (rho=64, n=50)
    # -------------------------------------------------------------
    ax_c = axes[1, 0]
    p_high = df_sum[(df_sum["category"] == "power") & (df_sum["rho"] == 64.0) & (df_sum["n_taxa"] == 50)].sort_values("theta")
    th_labels = [r"$\theta=5$", r"$\theta=10$", r"$\theta=20$"]
    x_th = np.arange(len(th_labels))

    for idx, (m_name, m_col, col) in enumerate(methods):
        rate_key = method_map[m_name]
        vals = p_high[rate_key].values
        ax_c.bar(x_th + (idx - 2) * width, vals, width, label=m_name, color=col, alpha=0.88)

    ax_c.set_xticks(x_th)
    ax_c.set_xticklabels(th_labels)
    ax_c.set_xlabel("Mutation Diversity Parameter")
    ax_c.set_ylabel("Sensitivity (%)")
    ax_c.set_title(r"C   Sensitivity at High Recombination ($\rho=64, n=50$)", loc="left", fontweight="bold")
    ax_c.set_ylim(0, 105)
    ax_c.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")

    # -------------------------------------------------------------
    # Panel D: Execution Latency (ms per alignment)
    # -------------------------------------------------------------
    ax_d = axes[1, 1]
    
    # Compare latency for n=10 vs n=50
    lat_data = [
        ("RhizAeon", df_raw[df_raw["n_taxa"] == 10]["rhiz_time_ms"].mean(), df_raw[df_raw["n_taxa"] == 50]["rhiz_time_ms"].mean(), COLORS["RhizAeon"]),
        ("3SEQ", df_raw[df_raw["n_taxa"] == 10]["three_time_ms"].mean(), df_raw[df_raw["n_taxa"] == 50]["three_time_ms"].mean(), COLORS["3SEQ"]),
        ("PHI Suite", df_raw[df_raw["n_taxa"] == 10]["phi_suite_time_ms"].mean(), df_raw[df_raw["n_taxa"] == 50]["phi_suite_time_ms"].mean(), COLORS["PHI"])
    ]

    x_tools = np.arange(len(lat_data))
    w_lat = 0.35
    ax_d.bar(x_tools - w_lat/2, [d[1] for d in lat_data], w_lat, label="n = 10 taxa", color="#94a3b8", edgecolor="none")
    ax_d.bar(x_tools + w_lat/2, [d[2] for d in lat_data], w_lat, label="n = 50 taxa", color="#0284c7", edgecolor="none")

    ax_d.set_xticks(x_tools)
    ax_d.set_xticklabels([d[0] for d in lat_data])
    ax_d.set_ylabel("Runtime Latency (ms)")
    ax_d.set_yscale("log")
    ax_d.set_title("D   Computational Scalability & Latency", loc="left", fontweight="bold")
    ax_d.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_d.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    # Convert PDF to high-resolution PNG using macOS sips
    try:
        import subprocess
        subprocess.run(["/usr/bin/sips", "-s", "format", "png", str(out_pdf), "--out", str(out_png)], check=True, stdout=subprocess.DEVNULL)
    except Exception:
        pass
    print(f"Saved publication figures to:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    sum_csv = base / "phipack_summary.csv"
    r_csv = base / "phipack_raw_results.csv"
    p_pdf = Path(__file__).resolve().parent.parent.parent / "paper" / "Figures" / "fig_phipack_replication.pdf"
    p_png = Path(__file__).resolve().parent.parent.parent / "paper" / "Figures" / "fig_phipack_replication.png"
    p_pdf.parent.mkdir(parents=True, exist_ok=True)
    if sum_csv.exists() and r_csv.exists():
        plot_phipack_figure(sum_csv, r_csv, p_pdf, p_png)
    else:
        print("Summary CSV not yet generated.")
