"""
benchmarks/multi_taxon_simulations/run_multi_taxon_baseline.py
==============================================================
Comprehensive Multi-Taxon Baseline Benchmark Harness for RhizAeon.

Executes all three multi-taxon options:
  Option A: Taxonomic Scaling (N in {5, 10, 20, 50, 100}, d in {2%, 5%})
  Option B: Ghost Donors / Incomplete Lineage Sampling (Delta d in {0%, 1%, 3%, 5%, 10%}, Complete Ghost)
  Option C: Circulating Recombinant Forms / Clade Expansion (M in {1, 2, 4, 8, 12} out of N=25)

Outputs trial-by-trial results, aggregated summary CSVs, and publication-grade figures.
"""

import os
import sys

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

import time
import argparse
from typing import Dict, Any, List, Tuple
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
import matplotlib.pyplot as plt

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, PROJECT_ROOT)

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
from benchmarks.multi_taxon_simulations.generator_multi import (
    generate_scaling_alignment, generate_ghost_alignment, generate_crf_alignment
)

OKABE_ITO = {
    "rhizaeon": "#D55E00",      # Vermillion
    "blue": "#0072B2",          # Sky Blue
    "green": "#009E73",         # Bluish Green
    "amber": "#E69F00",         # Amber
    "purple": "#CC79A7",        # Reddish Purple
    "gray": "#7F7F7F"
}

# -----------------------------------------------------------------------------
# Worker Functions for ProcessPoolExecutor
# -----------------------------------------------------------------------------

def evaluate_single_scaling_rep(args: Dict[str, Any]) -> Dict[str, Any]:
    N = args["N"]
    d = args["divergence"]
    rep_idx = args["rep_idx"]
    seed = args["seed"]
    is_null = args["is_null"]

    t0 = time.perf_counter()
    mat, taxa, L, meta = generate_scaling_alignment(
        N=N, divergence=d, tract_len=500, tract_start=1000, L=3000, seed=seed, is_null=is_null
    )
    t_gen = (time.perf_counter() - t0) * 1000

    # Scaled Z-score threshold for large N
    z_thresh = max(2.5, 2.0 + 0.25 * np.log(max(N, 3)))
    pir_thresh = 0.18

    t1 = time.perf_counter()
    eng = PrefixDistanceEngine(mat, codon_aligned=False)
    bps = run_recursive_partition_fda_screen(
        eng, taxa, min_flank=40, max_flank=150, step=25, min_z=z_thresh, min_pir=pir_thresh,
        frobenius_triage=True, max_depth=2, polish_ml=True
    )
    t_screen = (time.perf_counter() - t1) * 1000

    if is_null:
        return {
            "experiment": "scaling",
            "N": N,
            "divergence": d,
            "is_null": True,
            "rep_idx": rep_idx,
            "detected": len(bps) > 0,
            "n_bps": len(bps),
            "pwr_any": 0.0,
            "pwr_both": 0.0,
            "mae": np.nan,
            "runtime_ms": t_screen
        }
    else:
        true_bps = [1000, 1500]
        rec_bps = [b.breakpoint_nt for b in bps if b.recombinant_taxon == "Rec"]
        # Also check any detected bps
        all_bps = [b.breakpoint_nt for b in bps]
        target_bps = rec_bps if len(rec_bps) > 0 else all_bps

        hit1 = any(abs(b - true_bps[0]) <= 100 for b in target_bps)
        hit2 = any(abs(b - true_bps[1]) <= 100 for b in target_bps)

        pwr_any = float(hit1 or hit2)
        pwr_both = float(hit1 and hit2)

        maes = []
        for tb in true_bps:
            dists = [abs(b - tb) for b in target_bps if abs(b - tb) <= 150]
            if dists:
                maes.append(min(dists))
        mean_mae = float(np.mean(maes)) if maes else np.nan

        return {
            "experiment": "scaling",
            "N": N,
            "divergence": d,
            "is_null": False,
            "rep_idx": rep_idx,
            "detected": len(bps) > 0,
            "rec_matched": len(rec_bps) > 0,
            "n_bps": len(bps),
            "pwr_any": pwr_any,
            "pwr_both": pwr_both,
            "mae": mean_mae,
            "runtime_ms": t_screen
        }


def evaluate_single_ghost_rep(args: Dict[str, Any]) -> Dict[str, Any]:
    d_donor = args["donor_divergence"]
    is_ghost = args["is_complete_ghost"]
    rep_idx = args["rep_idx"]
    seed = args["seed"]

    mat, taxa, L, meta = generate_ghost_alignment(
        divergence=0.05, donor_divergence=d_donor, is_complete_ghost=is_ghost,
        tract_len=600, tract_start=1000, L=3000, seed=seed
    )

    t1 = time.perf_counter()
    eng = PrefixDistanceEngine(mat, codon_aligned=False)
    bps = run_recursive_partition_fda_screen(
        eng, taxa, min_flank=40, max_flank=150, step=25, min_z=2.5, min_pir=0.15,
        frobenius_triage=True, max_depth=2, polish_ml=True
    )
    t_screen = (time.perf_counter() - t1) * 1000

    true_bps = [1000, 1600]
    rec_bps = [b.breakpoint_nt for b in bps if b.recombinant_taxon == "Rec"]
    target_bps = rec_bps if rec_bps else [b.breakpoint_nt for b in bps]

    hit1 = any(abs(b - true_bps[0]) <= 100 for b in target_bps)
    hit2 = any(abs(b - true_bps[1]) <= 100 for b in target_bps)
    pwr_any = float(hit1 or hit2)
    pwr_both = float(hit1 and hit2)

    maes = []
    for tb in true_bps:
        dists = [abs(b - tb) for b in target_bps if abs(b - tb) <= 150]
        if dists:
            maes.append(min(dists))
    mean_mae = float(np.mean(maes)) if maes else np.nan

    top_z = max([b.kinetic_z for b in bps]) if bps else 0.0

    return {
        "experiment": "ghost",
        "donor_divergence": d_donor,
        "is_complete_ghost": is_ghost,
        "rep_idx": rep_idx,
        "detected": len(bps) > 0,
        "rec_matched": len(rec_bps) > 0,
        "pwr_any": pwr_any,
        "pwr_both": pwr_both,
        "mae": mean_mae,
        "max_z": top_z,
        "runtime_ms": t_screen
    }


def evaluate_single_crf_rep(args: Dict[str, Any]) -> Dict[str, Any]:
    M = args["M"]
    rep_idx = args["rep_idx"]
    seed = args["seed"]

    mat, taxa, L, meta = generate_crf_alignment(
        M_recombinants=M, N_total=25, divergence=0.05, clade_drift=0.005,
        tract_len=600, tract_start=1000, L=3000, seed=seed
    )

    t1 = time.perf_counter()
    eng = PrefixDistanceEngine(mat, codon_aligned=False)
    bps = run_recursive_partition_fda_screen(
        eng, taxa, min_flank=40, max_flank=150, step=25, min_z=1.8, min_pir=0.15,
        frobenius_triage=True, max_depth=2, polish_ml=True
    )
    t_screen = (time.perf_counter() - t1) * 1000

    true_bps = [1000, 1600]
    crf_names = set(meta["rec_taxa"])
    crf_bps = [b.breakpoint_nt for b in bps if b.recombinant_taxon in crf_names]
    crf_taxa_detected = set([b.recombinant_taxon for b in bps if b.recombinant_taxon in crf_names])

    target_bps = crf_bps if crf_bps else [b.breakpoint_nt for b in bps]

    hit1 = any(abs(b - true_bps[0]) <= 100 for b in target_bps)
    hit2 = any(abs(b - true_bps[1]) <= 100 for b in target_bps)
    pwr_any = float(hit1 or hit2)
    pwr_both = float(hit1 and hit2)

    clade_sensitivity = float(len(crf_taxa_detected) / M)

    maes = []
    for tb in true_bps:
        dists = [abs(b - tb) for b in target_bps if abs(b - tb) <= 150]
        if dists:
            maes.append(min(dists))
    mean_mae = float(np.mean(maes)) if maes else np.nan

    # Breakpoint concordance: std dev of predicted left and right breakpoints across clade members
    left_bps = [b.breakpoint_nt for b in bps if b.recombinant_taxon in crf_names and abs(b.breakpoint_nt - true_bps[0]) <= 150]
    std_left = float(np.std(left_bps)) if len(left_bps) > 1 else 0.0

    return {
        "experiment": "crf",
        "M": M,
        "rep_idx": rep_idx,
        "event_detected": pwr_any > 0,
        "clade_sensitivity": clade_sensitivity,
        "pwr_any": pwr_any,
        "pwr_both": pwr_both,
        "mae": mean_mae,
        "std_left_bp": std_left,
        "runtime_ms": t_screen
    }


# -----------------------------------------------------------------------------
# Main Execution Harness
# -----------------------------------------------------------------------------

def run_all_multi_taxon_benchmarks(n_workers: int = 3):
    print("================================================================================")
    print("MULTI-TAXON RHIZAEON BASELINE BENCHMARK")
    print(f"Evaluating Options A, B, and C with {n_workers} lightweight workers")
    print("================================================================================")

    # -------------------------------------------------------------------------
    # Option A: Taxonomic Scaling Tasks
    # -------------------------------------------------------------------------
    scaling_tasks = []
    for N in [5, 10, 20, 50, 100]:
        for d in [0.02, 0.05]:
            # Recombinant runs (15 reps)
            for r in range(15):
                scaling_tasks.append({
                    "N": N, "divergence": d, "rep_idx": r, "seed": 1000 + r * 11 + N * 7, "is_null": False
                })
            # Null runs (10 reps)
            for r in range(10):
                scaling_tasks.append({
                    "N": N, "divergence": d, "rep_idx": r, "seed": 2000 + r * 13 + N * 7, "is_null": True
                })

    print(f"Option A (Scaling): Launching {len(scaling_tasks)} tasks on {n_workers} workers...", flush=True)
    scaling_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(evaluate_single_scaling_rep, t) for t in scaling_tasks]
        for i, f in enumerate(as_completed(futures), 1):
            scaling_results.append(f.result())
            if i % 25 == 0 or i == len(scaling_tasks):
                print(f"  [Scaling] Progress: {i}/{len(scaling_tasks)} ({i/len(scaling_tasks)*100:.1f}%) completed", flush=True)

    df_scale = pd.DataFrame(scaling_results)
    df_scale.to_csv("benchmarks/multi_taxon_simulations/scaling_raw_results.csv", index=False)
    print("Option A Completed!\n", flush=True)

    # -------------------------------------------------------------------------
    # Option B: Ghost Donor Tasks
    # -------------------------------------------------------------------------
    ghost_tasks = []
    ghost_conditions = [
        (0.00, False),  # Exact donor
        (0.01, False),  # 1% donor drift
        (0.03, False),  # 3% donor drift
        (0.05, False),  # 5% donor drift
        (0.10, False),  # 10% donor drift
        (0.00, True)    # Complete Ghost (Clade B 100% missing)
    ]
    for d_donor, is_ghost in ghost_conditions:
        for r in range(20):
            ghost_tasks.append({
                "donor_divergence": d_donor, "is_complete_ghost": is_ghost,
                "rep_idx": r, "seed": 3000 + r * 17 + int(d_donor * 1000) + (100 if is_ghost else 0)
            })

    print(f"Option B (Ghost Donors): Launching {len(ghost_tasks)} tasks on {n_workers} workers...", flush=True)
    ghost_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(evaluate_single_ghost_rep, t) for t in ghost_tasks]
        for i, f in enumerate(as_completed(futures), 1):
            ghost_results.append(f.result())
            if i % 20 == 0 or i == len(ghost_tasks):
                print(f"  [Ghost] Progress: {i}/{len(ghost_tasks)} ({i/len(ghost_tasks)*100:.1f}%) completed", flush=True)

    df_ghost = pd.DataFrame(ghost_results)
    df_ghost.to_csv("benchmarks/multi_taxon_simulations/ghost_raw_results.csv", index=False)
    print("Option B Completed!\n", flush=True)

    # -------------------------------------------------------------------------
    # Option C: CRF Clade Expansion Tasks
    # -------------------------------------------------------------------------
    crf_tasks = []
    for M in [1, 2, 4, 8, 12]:
        for r in range(20):
            crf_tasks.append({
                "M": M, "rep_idx": r, "seed": 4000 + r * 19 + M * 23
            })

    print(f"Option C (CRFs): Launching {len(crf_tasks)} tasks on {n_workers} workers...", flush=True)
    crf_results = []
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = [executor.submit(evaluate_single_crf_rep, t) for t in crf_tasks]
        for i, f in enumerate(as_completed(futures), 1):
            crf_results.append(f.result())
            if i % 20 == 0 or i == len(crf_tasks):
                print(f"  [CRF] Progress: {i}/{len(crf_tasks)} ({i/len(crf_tasks)*100:.1f}%) completed", flush=True)

    df_crf = pd.DataFrame(crf_results)
    df_crf.to_csv("benchmarks/multi_taxon_simulations/crf_raw_results.csv", index=False)
    print("Option C Completed!\n", flush=True)

    # -------------------------------------------------------------------------
    # Summarize and Plot
    # -------------------------------------------------------------------------
    summarize_and_plot(df_scale, df_ghost, df_crf)


def summarize_and_plot(df_scale: pd.DataFrame, df_ghost: pd.DataFrame, df_crf: pd.DataFrame):
    # Summary A: Scaling
    rec_scale = df_scale[~df_scale["is_null"]]
    null_scale = df_scale[df_scale["is_null"]]

    sum_scale_rec = rec_scale.groupby(["N", "divergence"])[["pwr_any", "pwr_both", "mae", "runtime_ms"]].mean().reset_index()
    sum_scale_null = null_scale.groupby(["N", "divergence"])[["detected"]].mean().reset_index().rename(columns={"detected": "fpr"})
    sum_scale = pd.merge(sum_scale_rec, sum_scale_null, on=["N", "divergence"])
    sum_scale.to_csv("benchmarks/multi_taxon_simulations/multi_taxon_scaling_summary.csv", index=False)

    # Summary B: Ghost
    sum_ghost = df_ghost.groupby(["donor_divergence", "is_complete_ghost"])[["pwr_any", "pwr_both", "mae", "max_z"]].mean().reset_index()
    sum_ghost.to_csv("benchmarks/multi_taxon_simulations/multi_taxon_ghost_summary.csv", index=False)

    # Summary C: CRF
    sum_crf = df_crf.groupby(["M"])[["event_detected", "clade_sensitivity", "pwr_both", "mae", "std_left_bp"]].mean().reset_index()
    sum_crf.to_csv("benchmarks/multi_taxon_simulations/multi_taxon_crf_summary.csv", index=False)

    print("\n=== OPTION A: TAXONOMIC SCALING SUMMARY ===")
    print(sum_scale.to_string(index=False))

    print("\n=== OPTION B: GHOST DONOR SUMMARY ===")
    print(sum_ghost.to_string(index=False))

    print("\n=== OPTION C: CRF CLADE EXPANSION SUMMARY ===")
    print(sum_crf.to_string(index=False))

    # Generate Publication Figure
    plot_multi_taxon_figure(sum_scale, sum_ghost, sum_crf)


def plot_multi_taxon_figure(sum_scale: pd.DataFrame, sum_ghost: pd.DataFrame, sum_crf: pd.DataFrame):
    plt.rcParams["font.sans-serif"] = "Helvetica"
    plt.rcParams["axes.edgecolor"] = "#222222"
    plt.rcParams["axes.linewidth"] = 0.8

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), dpi=300)

    # Panel 1: Option A Scaling (Power & Runtime vs N)
    ax1 = axes[0]
    d05 = sum_scale[sum_scale["divergence"] == 0.05]
    d02 = sum_scale[sum_scale["divergence"] == 0.02]

    ax1.plot(d05["N"], d05["pwr_both"] * 100, "o-", color=OKABE_ITO["rhizaeon"], linewidth=2.0, label="Power (d = 5.0%)")
    ax1.plot(d02["N"], d02["pwr_both"] * 100, "s--", color=OKABE_ITO["amber"], linewidth=1.8, label="Power (d = 2.0%)")
    ax1.set_xlabel("Taxon Sample Size (N)", fontsize=9.5, fontweight="bold")
    ax1.set_ylabel("Dual-Breakpoint Power (%)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A. Taxonomic Scaling (N = 5 to 100)", loc="left", fontsize=10.5, fontweight="bold")
    ax1.set_ylim(-5, 105)
    ax1.legend(loc="lower left", fontsize=8.5)
    ax1.grid(alpha=0.25, linestyle=":")

    # Inset for Runtime
    ax1_in = ax1.inset_axes([0.55, 0.52, 0.40, 0.40])
    ax1_in.plot(d05["N"], d05["runtime_ms"], "d-", color=OKABE_ITO["blue"], linewidth=1.5)
    ax1_in.set_title("Runtime (ms)", fontsize=7.5)
    ax1_in.set_xlabel("N", fontsize=7)
    ax1_in.set_ylabel("ms", fontsize=7)
    ax1_in.tick_params(labelsize=6.5)
    ax1_in.grid(alpha=0.2, linestyle=":")

    # Panel 2: Option B Ghost Donors (Power vs Donor Drift)
    ax2 = axes[1]
    g_sampled = sum_ghost[~sum_ghost["is_complete_ghost"]].sort_values("donor_divergence")
    g_complete = sum_ghost[sum_ghost["is_complete_ghost"]]

    drift_pcts = g_sampled["donor_divergence"] * 100
    ax2.plot(drift_pcts, g_sampled["pwr_any"] * 100, "o-", color=OKABE_ITO["rhizaeon"], linewidth=2.0, label="Event Power")
    ax2.plot(drift_pcts, g_sampled["pwr_both"] * 100, "s--", color=OKABE_ITO["amber"], linewidth=1.8, label="Dual Breakpoint Power")

    # Add Complete Ghost point
    c_pwr_any = g_complete["pwr_any"].values[0] * 100
    c_pwr_both = g_complete["pwr_both"].values[0] * 100
    ax2.scatter([15], [c_pwr_any], color=OKABE_ITO["green"], s=80, marker="D", zorder=5, label=f"Complete Ghost ({c_pwr_any:.0f}%)")
    ax2.scatter([15], [c_pwr_both], color=OKABE_ITO["green"], s=80, marker="x", zorder=5)

    ax2.set_xlabel("Donor Lineage Divergence (%)", fontsize=9.5, fontweight="bold")
    ax2.set_ylabel("Detection Power (%)", fontsize=9.5, fontweight="bold")
    ax2.set_title("B. Ghost Donor & Missing Lineages", loc="left", fontsize=10.5, fontweight="bold")
    ax2.set_ylim(-5, 105)
    ax2.legend(loc="lower left", fontsize=8.5)
    ax2.grid(alpha=0.25, linestyle=":")

    # Panel 3: Option C CRF Clade Expansion (Clade Sensitivity & Concordance)
    ax3 = axes[2]
    ax3.plot(sum_crf["M"], sum_crf["event_detected"] * 100, "o-", color=OKABE_ITO["rhizaeon"], linewidth=2.0, label="Ancestral Event Recovery")
    ax3.plot(sum_crf["M"], sum_crf["clade_sensitivity"] * 100, "s--", color=OKABE_ITO["blue"], linewidth=1.8, label="Mean Clade Member Sensitivity")

    ax3.set_xlabel("CRF Clade Size (M Descendants out of N=25)", fontsize=9.5, fontweight="bold")
    ax3.set_ylabel("Sensitivity (%)", fontsize=9.5, fontweight="bold")
    ax3.set_title("C. CRF Clade Propagation", loc="left", fontsize=10.5, fontweight="bold")
    ax3.set_ylim(-5, 105)
    ax3.legend(loc="lower right", fontsize=8.5)
    ax3.grid(alpha=0.25, linestyle=":")

    plt.tight_layout()
    out_fig = "paper/figures/fig_multi_taxon_rhizaeon_baseline.png"
    plt.savefig(out_fig, dpi=300)
    plt.savefig(out_fig.replace(".png", ".pdf"))
    print(f"\nSaved publication figure to {out_fig} and .pdf")

if __name__ == "__main__":
    run_all_multi_taxon_benchmarks()
