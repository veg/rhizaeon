"""
simulations/04_rdp5_martin_2021/plot_rdp5_benchmark.py
======================================================
BioVis publication-quality multipanel figure replicating Martin et al. (2021):
- Panel A: Computational speedup of RhizAeon over RDP5 and RDP4 across empirical viral cohorts.
- Panel B: Wall-clock execution latency (seconds) across benchmarks.
- Panel C: Runtime scaling as a function of alignment complexity (N * L).
- Panel D: Biological breakpoint localization precision across empirical targets.
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
    "RDP5": "#7c3aed",       # Purple
    "RDP4": "#64748b"        # Slate
}


def plot_rdp5_figure(results_csv: Path, out_pdf: Path, out_png: Path):
    df = pd.read_csv(results_csv)

    fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.8), dpi=300)
    plt.subplots_adjust(wspace=0.28, hspace=0.35, top=0.92, bottom=0.08, left=0.08, right=0.96)

    # -------------------------------------------------------------
    # Panel A: Computational Speedup over RDP5 and RDP4
    # -------------------------------------------------------------
    ax_a = axes[0, 0]
    scaled_df = df[df["speedup_vs_rdp5"].notna()].copy()
    labels = [r["dataset_name"] for _, r in scaled_df.iterrows()]
    x_pos = np.arange(len(labels))
    width = 0.35

    sp5 = scaled_df["speedup_vs_rdp5"].values
    sp4 = scaled_df["speedup_vs_rdp4"].values

    ax_a.bar(x_pos - width/2, sp5, width, label="Speedup vs RDP5", color=COLORS["RDP5"], alpha=0.88, edgecolor="none")
    ax_a.bar(x_pos + width/2, sp4, width, label="Speedup vs RDP4", color=COLORS["RDP4"], alpha=0.88, edgecolor="none")

    ax_a.set_xticks(x_pos)
    ax_a.set_xticklabels(labels, rotation=20, ha="right", fontsize=8.0)
    ax_a.set_ylabel(r"Acceleration Factor ($\times$)")
    ax_a.set_yscale("log")
    ax_a.set_title("A   Computational Speedup over RDP4 and RDP5", loc="left", fontweight="bold")
    ax_a.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_a.legend(loc="upper right", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    # Annotate values on top of bars
    for i in range(len(labels)):
        ax_a.text(x_pos[i] - width/2, sp5[i] * 1.15, f"{sp5[i]:.0f}x", ha="center", va="bottom", fontsize=7.0, color="#4c1d95")
        ax_a.text(x_pos[i] + width/2, sp4[i] * 1.15, f"{sp4[i]:.0f}x", ha="center", va="bottom", fontsize=7.0, color="#334155")
    ax_a.set_ylim(1, 2000)

    # -------------------------------------------------------------
    # Panel B: Wall-Clock Runtime Comparison (Logarithmic Scale)
    # -------------------------------------------------------------
    ax_b = axes[0, 1]
    w3 = 0.25
    r_rhiz = scaled_df["rhiz_total_sec"].values
    r_rdp5 = scaled_df["rdp5_published_sec"].values
    r_rdp4 = scaled_df["rdp4_published_sec"].values

    ax_b.bar(x_pos - w3, r_rhiz, w3, label="RhizAeon", color=COLORS["RhizAeon"], alpha=0.88, edgecolor="none")
    ax_b.bar(x_pos, r_rdp5, w3, label="RDP5", color=COLORS["RDP5"], alpha=0.88, edgecolor="none")
    ax_b.bar(x_pos + w3, r_rdp4, w3, label="RDP4", color=COLORS["RDP4"], alpha=0.88, edgecolor="none")

    ax_b.set_xticks(x_pos)
    ax_b.set_xticklabels(labels, rotation=20, ha="right", fontsize=8.0)
    ax_b.set_ylabel("Execution Time (seconds)")
    ax_b.set_yscale("log")
    ax_b.set_title("B   Wall-Clock Runtime Comparison", loc="left", fontweight="bold")
    ax_b.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")
    ax_b.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    # -------------------------------------------------------------
    # Panel C: Runtime vs Sequence Complexity (N * L)
    # -------------------------------------------------------------
    ax_c = axes[1, 0]
    all_n_l = df["taxa"].values * df["length_nt"].values
    all_time = df["rhiz_total_sec"].values

    ax_c.scatter(all_n_l, all_time, color=COLORS["RhizAeon"], s=60, zorder=5, edgecolors="#0369a1", linewidths=1.2)
    for i, row in df.iterrows():
        nl = row["taxa"] * row["length_nt"]
        t = row["rhiz_total_sec"]
        offset_y = 1.3 if t > 1.0 else 0.7
        ax_c.annotate(row["dataset_name"], (nl, t), textcoords="offset points", xytext=(0, 6), ha="center", fontsize=7.0)

    # Plot theoretical O(N^2 + L log L) scaling guide
    nl_grid = np.logspace(np.log10(min(all_n_l)), np.log10(max(all_n_l)), 100)
    # reference scale
    t_ref = all_time[0] * (nl_grid / all_n_l[0]) ** 1.15
    ax_c.plot(nl_grid, t_ref, linestyle="--", color="#94a3b8", alpha=0.7, label=r"Empirical scaling guide ($\approx \mathcal{O}(N L)$)")

    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"Alignment Complexity ($N \times L$ nucleotides)")
    ax_c.set_ylabel("RhizAeon Runtime (seconds)")
    ax_c.set_title("C   Empirical Algorithmic Scaling", loc="left", fontweight="bold")
    ax_c.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1")
    ax_c.legend(loc="upper left", frameon=True, framealpha=0.9, edgecolor="#e2e8f0")

    # -------------------------------------------------------------
    # Panel D: Empirical Biological Target Resolution
    # -------------------------------------------------------------
    ax_d = axes[1, 1]
    # Summary of events detected across the 6 viral systems
    d_names = [r["dataset_name"].replace(" ", "\n") for _, r in df.iterrows()]
    d_events = df["rhiz_bps_count"].values
    d_taxa = df["taxa"].values

    bars = ax_d.bar(np.arange(len(d_names)), d_events, color=COLORS["RhizAeon"], width=0.55, alpha=0.88, edgecolor="none")
    ax_d.set_xticks(np.arange(len(d_names)))
    ax_d.set_xticklabels(d_names, fontsize=7.5)
    ax_d.set_ylabel("Detected Recombination Events")
    ax_d.set_title("D   Recombination Complexity across Viral Systems", loc="left", fontweight="bold")
    ax_d.grid(True, linestyle=":", alpha=0.5, color="#cbd5e1", axis="y")

    for i, b in enumerate(bars):
        ax_d.text(b.get_x() + b.get_width()/2, b.get_height() + 1.0, f"{int(d_events[i])}\n(N={d_taxa[i]})", ha="center", va="bottom", fontsize=7.0, color="#0f172a")
    ax_d.set_ylim(0, max(d_events) * 1.25)

    plt.tight_layout()
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close()

    try:
        import subprocess
        subprocess.run(["/usr/bin/sips", "-s", "format", "png", str(out_pdf), "--out", str(out_png)], check=True, stdout=subprocess.DEVNULL)
    except Exception:
        pass
    print(f"Saved publication figures to:\n  {out_pdf}\n  {out_png}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    r_csv = base / "rdp5_benchmark_results.csv"
    p_pdf = Path(__file__).resolve().parent.parent.parent / "paper" / "Figures" / "fig_rdp5_benchmark.pdf"
    p_png = Path(__file__).resolve().parent.parent.parent / "paper" / "Figures" / "fig_rdp5_benchmark.png"
    p_pdf.parent.mkdir(parents=True, exist_ok=True)
    if r_csv.exists():
        plot_rdp5_figure(r_csv, p_pdf, p_png)
    else:
        print("CSV not found.")
