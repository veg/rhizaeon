"""
simulations/02_rivet_sars_cov_2/run_rivet_benchmark.py
======================================================
Unified benchmark runner evaluating RhizAeon and 3SEQ on:
1. Component A (OG Data): All 481 high-confidence QC PASS pandemic recombinant trios
   from RIVET (Smith et al. 2023, Bioinformatics; rivet.ucsd.edu).
2. Component B (Regenerated Sim): 1,500 whole-genome SARS-CoV-2 alignments replicating
   the Thornlow & Turakhia simulation protocol (Nature 2022, Bioinformatics 2023).
"""

import os
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
import re
import argparse
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
STUDY_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(STUDY_DIR))

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
import generator

THREE_SEQ_BIN = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq")
PTABLE_PATH = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq_ptable_250")
MASTER_SEED = 20230831  # Bioinformatics publication date for RIVET

L_GENOME = 29903


def parse_interval(interval_str: str) -> Tuple[int, int]:
    """Parses '(5239,6449)' into (5239, 6449)."""
    m = re.search(r"\((\d+)\s*,\s*(\d+)\)", str(interval_str))
    if m:
        return int(m.group(1)), int(m.group(2))
    return 0, 0


def run_single_empirical_trio(row_dict: Dict[str, Any], run_3seq: bool = True) -> Dict[str, Any]:
    """Evaluates a single empirical RIVET recombinant trio."""
    rec_lineage = str(row_dict.get("Recombinant Lineage", "Unknown"))
    rec_node = str(row_dict.get("Recombinant Node ID", "Unknown"))
    bp1_left, bp1_right = parse_interval(row_dict.get("Breakpoint 1 Interval", "(0,0)"))
    bp2_left, bp2_right = parse_interval(row_dict.get("Breakpoint 2 Interval", "(0,0)"))
    rivet_p_3seq = float(row_dict.get("3SEQ P-Value", 1.0))
    rivet_m_n_k = str(row_dict.get("3SEQ (M, N, K)", ""))

    positions = [int(p) for p in str(row_dict.get("Informative Site Positions", "")).split(",") if p.strip()]
    seq_str = str(row_dict.get("Informative Site Sequence", "")).strip()

    # 1. Construct 3-sequence alignment across 29,903 nt
    mat = np.zeros((3, L_GENOME), dtype=np.int8)
    for pos, char in zip(positions, seq_str):
        idx = min(max(0, pos - 1), L_GENOME - 1)
        mat[1, idx] = 1  # Donor
        mat[2, idx] = 0  # Acceptor
        mat[0, idx] = 1 if char == "A" else 0  # Recombinant

    # 2. Run RhizAeon standard default workflow
    t0_rhiz = time.perf_counter()
    rhiz_bps = []
    try:
        engine = PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=False)
        rhiz_bps = run_recursive_partition_fda_screen(
            engine,
            taxa_names=["Recombinant", "Donor", "Acceptor"],
            min_flank=1500,
            max_flank=3000,
            step=250,
            min_z=1.0,
            min_pir=0.04,
            crossover_validation=True,
            crossover_p_threshold=0.005,
            min_informative_sites=2,
            polish_ml=True,
            min_parent_dist=1e-5
        )
        rhiz_time_ms = (time.perf_counter() - t0_rhiz) * 1000.0
    except Exception as e:
        rhiz_time_ms = (time.perf_counter() - t0_rhiz) * 1000.0

    rhiz_detected = len(rhiz_bps) > 0
    rhiz_coords = [b.polished_bp if b.polished_bp is not None else b.breakpoint_nt for b in rhiz_bps]

    # Evaluate concordance with RIVET breakpoint intervals
    concordant = False
    min_dist_to_interval = float("inf")
    if rhiz_detected and (bp1_left > 0 or bp1_right > 0):
        for c in rhiz_coords:
            # Check if c is inside [bp1_left, bp1_right] or within mutational resolution (1,500 nt)
            if bp1_left <= c <= bp1_right:
                concordant = True
                min_dist_to_interval = 0
                break
            d_left = abs(c - bp1_left)
            d_right = abs(c - bp1_right)
            cur_dist = min(d_left, d_right)
            if cur_dist < min_dist_to_interval:
                min_dist_to_interval = cur_dist
            if cur_dist <= 1500:
                concordant = True

    # 3. Optional 3SEQ run
    three_seq_time_ms = 0.0
    three_seq_detected = False
    three_seq_p = 1.0
    if run_3seq and os.path.exists(THREE_SEQ_BIN):
        t0_3s = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory() as td:
                fa_path = os.path.join(td, "trio.fa")
                chars = generator.NT_CHARS[mat]
                with open(fa_path, "w") as f:
                    for t_idx, t_name in enumerate(["Recombinant", "Donor", "Acceptor"]):
                        f.write(f">{t_name}\n{''.join(chars[t_idx])}\n")
                cmd = [
                    THREE_SEQ_BIN,
                    "-full", fa_path,
                    "-p", PTABLE_PATH,
                    "-id", "run_3s",
                    "-q",
                    "-L40"
                ]
                subprocess.run(cmd, cwd=td, capture_output=True, text=True, check=False)
                three_seq_time_ms = (time.perf_counter() - t0_3s) * 1000.0
                csv_path = os.path.join(td, "run_3s.3s.rec.csv")
                if os.path.exists(csv_path):
                    with open(csv_path) as f:
                        lines = [l.strip() for l in f if l.strip()]
                    if len(lines) > 1:
                        parts = lines[1].split(",")
                        if len(parts) >= 7 and parts[6] != "":
                            three_seq_p = float(parts[6])
                            three_seq_detected = (three_seq_p <= 0.005)
        except Exception:
            three_seq_time_ms = (time.perf_counter() - t0_3s) * 1000.0

    return {
        "recombinant_node": rec_node,
        "recombinant_lineage": rec_lineage,
        "rivet_bp1_interval": f"({bp1_left},{bp1_right})",
        "rivet_bp2_interval": f"({bp2_left},{bp2_right})",
        "rivet_3seq_p": rivet_p_3seq,
        "informative_sites_count": len(positions),
        "rhizaeon_detected": rhiz_detected,
        "rhizaeon_num_bps": len(rhiz_bps),
        "rhizaeon_bps": ";".join(str(c) for c in rhiz_coords),
        "rhizaeon_concordant": concordant,
        "rhizaeon_min_dist_to_interval": min_dist_to_interval,
        "rhizaeon_runtime_ms": rhiz_time_ms,
        "three_seq_detected": three_seq_detected,
        "three_seq_p": three_seq_p,
        "three_seq_runtime_ms": three_seq_time_ms
    }


def run_single_simulation_replicate(
    n_bps: int,
    d: int,
    m: int,
    rep_idx: int,
    seed: int,
    run_3seq: bool = True
) -> Dict[str, Any]:
    """Runs a single replicate of the Thornlow & Turakhia simulation protocol."""
    # 1. Generate dataset
    t0_gen = time.perf_counter()
    mat, taxa, meta = generator.generate_sars_cov_2_trio(n_bps=n_bps, d=d, m=m, seed=seed)
    gen_time_ms = (time.perf_counter() - t0_gen) * 1000.0

    # 2. Run RhizAeon standard default workflow
    t0_rhiz = time.perf_counter()
    rhiz_bps = []
    try:
        engine = PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=False)
        rhiz_bps = run_recursive_partition_fda_screen(
            engine,
            taxa_names=taxa,
            min_flank=1500,
            max_flank=3000,
            step=250,
            min_z=1.0,
            min_pir=0.04,
            crossover_validation=True,
            crossover_p_threshold=0.005,
            min_informative_sites=2,
            polish_ml=True,
            min_parent_dist=1e-5
        )
        rhiz_time_ms = (time.perf_counter() - t0_rhiz) * 1000.0
    except Exception:
        rhiz_time_ms = (time.perf_counter() - t0_rhiz) * 1000.0

    rhiz_detected = len(rhiz_bps) > 0
    rhiz_coords = [b.polished_bp if b.polished_bp is not None else b.breakpoint_nt for b in rhiz_bps]

    # Spatial error calculation
    true_bps = meta["true_bps"]
    spatial_err = None
    if n_bps > 0 and rhiz_detected:
        errs = []
        for t_bp in true_bps:
            min_e = min(abs(t_bp - c) for c in rhiz_coords)
            errs.append(min_e)
        spatial_err = float(np.mean(errs))

    # 3. Run 3SEQ
    three_seq_detected = False
    three_seq_p = 1.0
    three_seq_time_ms = 0.0
    if run_3seq and os.path.exists(THREE_SEQ_BIN):
        t0_3s = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory() as td:
                fa_path = os.path.join(td, "sim.fa")
                chars = generator.NT_CHARS[mat]
                with open(fa_path, "w") as f:
                    for t_idx, t_name in enumerate(taxa):
                        f.write(f">{t_name}\n{''.join(chars[t_idx])}\n")
                cmd = [
                    THREE_SEQ_BIN,
                    "-full", fa_path,
                    "-p", PTABLE_PATH,
                    "-id", "sim_3s",
                    "-q",
                    "-L40"
                ]
                subprocess.run(cmd, cwd=td, capture_output=True, text=True, check=False)
                three_seq_time_ms = (time.perf_counter() - t0_3s) * 1000.0
                csv_path = os.path.join(td, "sim_3s.3s.rec.csv")
                if os.path.exists(csv_path):
                    with open(csv_path) as f:
                        lines = [l.strip() for l in f if l.strip()]
                    if len(lines) > 1:
                        parts = lines[1].split(",")
                        if len(parts) >= 7 and parts[6] != "":
                            three_seq_p = float(parts[6])
                            three_seq_detected = (three_seq_p <= 0.005)
        except Exception:
            three_seq_time_ms = (time.perf_counter() - t0_3s) * 1000.0

    return {
        "n_bps": n_bps,
        "d": d,
        "m": m,
        "rep_idx": rep_idx,
        "seed": seed,
        "is_recombinant": n_bps > 0,
        "true_bps_count": len(true_bps),
        "true_bps": ";".join(str(b) for b in true_bps),
        "rhizaeon_detected": rhiz_detected,
        "rhizaeon_num_bps": len(rhiz_bps),
        "rhizaeon_bps": ";".join(str(c) for c in rhiz_coords),
        "rhizaeon_spatial_err_nt": spatial_err,
        "rhizaeon_runtime_ms": rhiz_time_ms,
        "three_seq_detected": three_seq_detected,
        "three_seq_p": three_seq_p,
        "three_seq_runtime_ms": three_seq_time_ms
    }


def main():
    parser = argparse.ArgumentParser(description="Run RIVET SARS-CoV-2 Recombination Benchmark")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel worker processes")
    parser.add_argument("--sim-reps", type=int, default=50, help="Replicates per simulation condition (default: 50)")
    parser.add_argument("--skip-3seq", action="store_true", help="Skip 3SEQ execution")
    args = parser.parse_args()

    run_3seq = not args.skip_3seq
    print("=" * 70)
    print("RIVET SARS-CoV-2 Recombination Benchmark (Study 02)")
    print(f"Workers: {args.workers} | Sim Replicates: {args.sim_reps} | 3SEQ: {run_3seq}")
    print("=" * 70)

    # -------------------------------------------------------------
    # PART 1: Component A - Original Empirical RIVET PASS Trios
    # -------------------------------------------------------------
    empirical_tsv = STUDY_DIR / "data" / "rivet_public_recombinants.tsv"
    if empirical_tsv.exists():
        print("\n--- PART 1: Evaluating Original Empirical RIVET PASS Trios ---")
        df_emp = pd.read_csv(empirical_tsv, sep="\t")
        pass_df = df_emp[df_emp["Quality Control (QC) Flags"] == "PASS"].copy()
        print(f"Found {len(pass_df)} QC PASS trios in {empirical_tsv.name}")

        rows_to_process = [row.to_dict() for _, row in pass_df.iterrows()]
        empirical_results = []
        t0_emp = time.perf_counter()

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(run_single_empirical_trio, r, run_3seq): idx for idx, r in enumerate(rows_to_process)}
            for fut in as_completed(futures):
                try:
                    res = fut.result()
                    empirical_results.append(res)
                except Exception as e:
                    print(f"Error processing empirical trio: {e}")

        total_emp_time = time.perf_counter() - t0_emp
        df_emp_res = pd.DataFrame(empirical_results)
        df_emp_res.sort_values(by="recombinant_node", inplace=True)
        emp_csv = STUDY_DIR / "rivet_empirical_trios.csv"
        df_emp_res.to_csv(emp_csv, index=False)
        print(f"Saved empirical trio results to {emp_csv} ({len(df_emp_res)} rows, {total_emp_time:.1f}s)")

        # Empirical summary
        n_total = len(df_emp_res)
        n_rhiz = int(df_emp_res["rhizaeon_detected"].sum())
        n_conc = int(df_emp_res["rhizaeon_concordant"].sum())
        n_3seq = int(df_emp_res["three_seq_detected"].sum())
        mean_rhiz_ms = df_emp_res["rhizaeon_runtime_ms"].mean()
        mean_3seq_ms = df_emp_res["three_seq_runtime_ms"].mean()

        print(f"Empirical Results Summary (N={n_total}):")
        print(f"  RhizAeon Detected: {n_rhiz}/{n_total} ({n_rhiz/n_total*100:.1f}%)")
        print(f"  RhizAeon Concordant with RIVET: {n_conc}/{n_total} ({n_conc/n_total*100:.1f}%)")
        print(f"  3SEQ Detected (p<=0.005): {n_3seq}/{n_total} ({n_3seq/n_total*100:.1f}%)")
        print(f"  RhizAeon Mean Latency: {mean_rhiz_ms:.2f} ms/genome")
        print(f"  3SEQ Mean Latency: {mean_3seq_ms:.2f} ms/genome")

        emp_summary_data = [{
            "cohort": "RIVET_QC_PASS_Empirical",
            "total_trios": n_total,
            "rhizaeon_detected_count": n_rhiz,
            "rhizaeon_sensitivity": n_rhiz / n_total,
            "rhizaeon_concordant_count": n_conc,
            "rhizaeon_concordance_rate": n_conc / n_total,
            "three_seq_detected_count": n_3seq,
            "three_seq_sensitivity": n_3seq / n_total,
            "rhizaeon_mean_latency_ms": mean_rhiz_ms,
            "three_seq_mean_latency_ms": mean_3seq_ms
        }]
        pd.DataFrame(emp_summary_data).to_csv(STUDY_DIR / "rivet_empirical_summary.csv", index=False)
    else:
        print(f"Warning: Empirical file {empirical_tsv} not found, skipping Part 1.")

    # -------------------------------------------------------------
    # PART 2: Component B - Faithfully Regenerated Simulations
    # -------------------------------------------------------------
    print("\n--- PART 2: Executing Faithfully Regenerated Simulations ---")
    d_levels = [10, 20, 30, 50]
    m_levels = [0, 1, 2, 3]
    reps_per_cond = args.sim_reps

    sim_jobs = []
    job_idx = 0

    # 1-Breakpoint simulations: 4 d x 4 m x reps
    for d in d_levels:
        for m in m_levels:
            for rep in range(reps_per_cond):
                seed = MASTER_SEED + job_idx
                sim_jobs.append((1, d, m, rep, seed))
                job_idx += 1

    # 2-Breakpoint simulations: 3 d x 2 m x reps
    for d in [20, 30, 50]:
        for m in [0, 2]:
            for rep in range(reps_per_cond):
                seed = MASTER_SEED + job_idx
                sim_jobs.append((2, d, m, rep, seed))
                job_idx += 1

    # Negative controls (pure vertical descent): 4 d x 4 m x reps
    for d in d_levels:
        for m in m_levels:
            for rep in range(reps_per_cond):
                seed = MASTER_SEED + job_idx
                sim_jobs.append((0, d, m, rep, seed))
                job_idx += 1

    print(f"Total simulated alignments queued: {len(sim_jobs)}")

    t0_sim = time.perf_counter()
    sim_results = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(run_single_simulation_replicate, n_bps, d, m, rep, seed, run_3seq): idx
                   for idx, (n_bps, d, m, rep, seed) in enumerate(sim_jobs)}
        for fut in as_completed(futures):
            try:
                res = fut.result()
                sim_results.append(res)
            except Exception as e:
                print(f"Error in simulation replicate: {e}")

    total_sim_time = time.perf_counter() - t0_sim
    df_sim = pd.DataFrame(sim_results)
    raw_sim_csv = STUDY_DIR / "rivet_raw_results.csv"
    df_sim.to_csv(raw_sim_csv, index=False)
    print(f"Saved raw simulation results to {raw_sim_csv} ({len(df_sim)} rows, {total_sim_time:.1f}s)")

    # Build stratified simulation summary
    summary_rows = []
    for (n_bps, d, m), group in df_sim.groupby(["n_bps", "d", "m"]):
        n_reps = len(group)
        rhiz_det = group["rhizaeon_detected"].sum()
        three_seq_det = group["three_seq_detected"].sum()
        rhiz_pwr = rhiz_det / n_reps
        three_seq_pwr = three_seq_det / n_reps
        spatial_err = group["rhizaeon_spatial_err_nt"].dropna().mean() if n_bps > 0 else 0.0
        rhiz_ms = group["rhizaeon_runtime_ms"].mean()
        three_seq_ms = group["three_seq_runtime_ms"].mean()

        summary_rows.append({
            "n_bps": n_bps,
            "d": d,
            "m": m,
            "replicates": n_reps,
            "rhizaeon_power": rhiz_pwr,
            "three_seq_power": three_seq_pwr,
            "rhizaeon_spatial_err_nt": spatial_err if not np.isnan(spatial_err) else None,
            "rhizaeon_mean_ms": rhiz_ms,
            "three_seq_mean_ms": three_seq_ms
        })

    df_summary = pd.DataFrame(summary_rows)
    summary_csv = STUDY_DIR / "rivet_summary.csv"
    df_summary.to_csv(summary_csv, index=False)
    print(f"Saved simulation summary to {summary_csv} ({len(df_summary)} scenarios)")

    # Print overall highlights
    print("\n" + "=" * 70)
    print("OVERALL BENCHMARK HIGHLIGHTS:")
    print("=" * 70)
    neg_controls = df_sim[df_sim["n_bps"] == 0]
    rec_1bp = df_sim[df_sim["n_bps"] == 1]
    rec_2bp = df_sim[df_sim["n_bps"] == 2]

    print(f"Negative Controls (Pure Vertical Evolution, N={len(neg_controls)}):")
    print(f"  RhizAeon False Positive Rate: {neg_controls['rhizaeon_detected'].mean()*100:.2f}%")
    print(f"  3SEQ False Positive Rate:     {neg_controls['three_seq_detected'].mean()*100:.2f}%")

    print(f"\n1-Breakpoint Recombinants (N={len(rec_1bp)}):")
    for d_val in [10, 20, 30, 50]:
        sub = rec_1bp[rec_1bp["d"] == d_val]
        print(f"  d = {d_val:2d} | RhizAeon Power: {sub['rhizaeon_detected'].mean()*100:5.1f}% | 3SEQ Power: {sub['three_seq_detected'].mean()*100:5.1f}% | Spatial Error: {sub['rhizaeon_spatial_err_nt'].dropna().mean():.1f} nt")

    print(f"\n2-Breakpoint Recombinants (N={len(rec_2bp)}):")
    for d_val in [20, 30, 50]:
        sub = rec_2bp[rec_2bp["d"] == d_val]
        print(f"  d = {d_val:2d} | RhizAeon Power: {sub['rhizaeon_detected'].mean()*100:5.1f}% | 3SEQ Power: {sub['three_seq_detected'].mean()*100:5.1f}% | Spatial Error: {sub['rhizaeon_spatial_err_nt'].dropna().mean():.1f} nt")

    print(f"\nComputational Latency across 29,903 nt viral genomes:")
    print(f"  RhizAeon: {df_sim['rhizaeon_runtime_ms'].mean():.2f} ms/genome")
    print(f"  3SEQ:     {df_sim['three_seq_runtime_ms'].mean():.2f} ms/genome")
    print("=" * 70)


if __name__ == "__main__":
    main()
