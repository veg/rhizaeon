"""
benchmarks/posada_2001_benchmark/rerun_posada_modern_rhizaeon.py
================================================================
Re-runs the full 3,200 dataset Posada & Crandall (2001) benchmark suite
using the latest calibrated RhizAeon framework (calibrated crossover validation,
profile ML polisher, and prefix distance tensor engine).
"""

import os
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
import math
from pathlib import Path
from typing import Dict, Any, Tuple, List
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
from benchmarks.posada_2001_benchmark.generator import generate_posada_dataset
from benchmarks.posada_2001_benchmark.scenarios import build_posada_scenarios


def evaluate_single_posada_run(args: Tuple[int, str, str, float, float, Any, int, bool, int, int, bool, float]) -> Dict[str, Any]:
    """Evaluates a single recorded Posada replicate using modern RhizAeon."""
    row_idx, sc_id, cat, theta, rho, alpha_val, seed, true_is_rec, true_parts, true_n_bps, three_det, three_time = args
    alpha = None if str(alpha_val) in ("inf", "None") else float(alpha_val)

    # 1. Re-generate identical dataset using recorded seed
    mat, taxa, meta = generate_posada_dataset(
        n=10, L=1000, theta=theta, rho=rho, alpha=alpha, seed=seed
    )
    true_bps = meta.get("true_breakpoints", [])

    # 2. Execute modern RhizAeon
    t0_rhiz = time.perf_counter()
    try:
        engine = PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=True)
        bps = run_recursive_partition_fda_screen(
            engine,
            taxa_names=taxa,
            min_z=1.8,
            min_pir=0.06,
            crossover_validation=True,
            crossover_p_threshold=0.005,
            min_informative_sites=3,
            polish_ml=True
        )
        rhiz_time_ms = (time.perf_counter() - t0_rhiz) * 1000.0
        rhiz_detected = len(bps) > 0
        rhiz_num_bps = len(bps)
        rhiz_bps_coords = [b.polished_bp if b.polished_bp is not None else b.breakpoint_nt for b in bps]
        rhiz_bps_str = ";".join(str(c) for c in rhiz_bps_coords)
    except Exception as e:
        rhiz_time_ms = 0.0
        rhiz_detected = False
        rhiz_num_bps = 0
        rhiz_bps_coords = []
        rhiz_bps_str = f"ERROR: {e}"

    # Spatial breakpoint accuracy
    spatial_errs = []
    if true_is_rec and len(true_bps) > 0 and len(rhiz_bps_coords) > 0:
        for t_bp in true_bps:
            min_e = min(abs(p_bp - t_bp) for p_bp in rhiz_bps_coords)
            spatial_errs.append(min_e)
    mean_spatial_err = float(np.mean(spatial_errs)) if spatial_errs else -1.0

    return {
        "row_idx": row_idx,
        "scenario_id": sc_id,
        "category": cat,
        "theta": theta,
        "rho": rho,
        "alpha": alpha_val,
        "seed": seed,
        "true_is_recombinant": true_is_rec,
        "true_num_partitions": true_parts,
        "true_num_breakpoints": true_n_bps,
        "rhizaeon_detected": rhiz_detected,
        "rhizaeon_num_bps": rhiz_num_bps,
        "rhizaeon_bps": rhiz_bps_str,
        "rhizaeon_spatial_err_nt": mean_spatial_err,
        "rhizaeon_runtime_ms": round(rhiz_time_ms, 2),
        "three_seq_detected": three_det,
        "three_seq_runtime_ms": three_time
    }


def execute_posada_modern_rerun(max_workers: int = 10):
    raw_csv = Path("benchmarks/posada_2001_benchmark/posada_raw_results.csv")
    df = pd.read_csv(raw_csv)
    total = len(df)
    print(f"[*] Re-evaluating all {total} Posada & Crandall datasets using latest RhizAeon code across {max_workers} workers...")

    task_args = []
    for idx, row in df.iterrows():
        task_args.append((
            idx,
            row["scenario_id"],
            row["category"],
            float(row["theta"]),
            float(row["rho"]),
            row["alpha"],
            int(row["seed"]),
            bool(row["true_is_recombinant"]),
            int(row["true_num_partitions"]),
            int(row["true_num_breakpoints"]),
            bool(row["three_seq_detected"]),
            float(row["three_seq_runtime_ms"])
        ))

    t0 = time.perf_counter()
    results = [None] * total
    completed = 0

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(evaluate_single_posada_run, arg): arg[0] for arg in task_args}
        for fut in as_completed(futures):
            res = fut.result()
            results[res["row_idx"]] = res
            completed += 1
            if completed % 400 == 0 or completed == total:
                elapsed = time.perf_counter() - t0
                print(f"[{completed}/{total}] ({completed/total*100:.1f}%) Completed in {elapsed:.1f}s")

    # Update raw DataFrame
    new_df = pd.DataFrame(results)
    new_df.to_csv(raw_csv, index=False)
    print(f"[*] Saved updated raw results to {raw_csv}")

    # Generate updated summary table
    summary_records = []
    for sc_id, group in new_df.groupby("scenario_id", sort=False):
        first = group.iloc[0]
        n_reps = len(group)
        rhiz_det = float(group["rhizaeon_detected"].mean() * 100.0)
        three_det = float(group["three_seq_detected"].mean() * 100.0)
        rhiz_time = float(group["rhizaeon_runtime_ms"].mean())
        three_time = float(group["three_seq_runtime_ms"].mean())
        mean_parts = float(group["true_num_partitions"].mean())
        
        valid_errs = group[group["rhizaeon_spatial_err_nt"] >= 0]["rhizaeon_spatial_err_nt"]
        mean_spat_err = float(valid_errs.mean()) if len(valid_errs) > 0 else -1.0

        summary_records.append({
            "scenario_id": sc_id,
            "category": first["category"],
            "theta": first["theta"],
            "rho": first["rho"],
            "alpha": first["alpha"],
            "n_reps": n_reps,
            "mean_true_partitions": mean_parts,
            "rhizaeon_detection_pct": rhiz_det,
            "three_seq_detection_pct": three_det,
            "rhizaeon_mean_spatial_err_nt": mean_spat_err,
            "rhizaeon_mean_time_ms": round(rhiz_time, 2),
            "three_seq_mean_time_ms": round(three_time, 2)
        })

    summary_df = pd.DataFrame(summary_records)
    summary_csv = Path("benchmarks/posada_2001_benchmark/posada_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    print(f"[*] Saved updated summary table to {summary_csv}")


if __name__ == "__main__":
    execute_posada_modern_rerun(max_workers=10)
