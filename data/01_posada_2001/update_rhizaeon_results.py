"""
benchmarks/posada_2001_benchmark/update_rhizaeon_results.py
===========================================================
Parallel evaluator that re-evaluates RhizAeon across the 3,200 simulated datasets
from posada_raw_results.csv using the exact recorded seeds, and updates
posada_raw_results.csv and posada_summary.csv.
"""

import os
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
from pathlib import Path
from typing import Dict, Any, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
from benchmarks.posada_2001_benchmark.generator import generate_posada_dataset
from benchmarks.posada_2001_benchmark.scenarios import build_posada_scenarios


def eval_rhizaeon_single(args: Tuple[int, float, float, Any, int]) -> Dict[str, Any]:
    row_idx, theta, rho, alpha_val, seed = args
    alpha = None if str(alpha_val) == "inf" else float(alpha_val)
    
    mat, taxa, meta = generate_posada_dataset(
        n=10, L=1000, theta=theta, rho=rho, alpha=alpha, seed=seed
    )
    
    t0 = time.perf_counter()
    try:
        engine = PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=True)
        bps = run_recursive_partition_fda_screen(
            engine,
            taxa_names=taxa,
            min_z=1.8,
            min_pir=0.08,
            crossover_validation=True,
            crossover_p_threshold=0.05
        )
        t_ms = (time.perf_counter() - t0) * 1000.0
        return {
            "row_idx": row_idx,
            "rhizaeon_detected": len(bps) > 0,
            "rhizaeon_num_bps": len(bps),
            "rhizaeon_bps": ";".join(str(b.polished_bp) for b in bps),
            "rhizaeon_runtime_ms": round(t_ms, 2)
        }
    except Exception as e:
        return {
            "row_idx": row_idx,
            "rhizaeon_detected": False,
            "rhizaeon_num_bps": 0,
            "rhizaeon_bps": f"ERROR: {e}",
            "rhizaeon_runtime_ms": 0.0
        }


def main():
    csv_path = Path("benchmarks/posada_2001_benchmark/posada_raw_results.csv")
    df = pd.read_csv(csv_path)
    total = len(df)
    print(f"[*] Dispatching parallel RhizAeon evaluation for {total} datasets across 3 workers...")

    tasks = [
        (idx, float(row["theta"]), float(row["rho"]), row["alpha"], int(row["seed"]))
        for idx, row in df.iterrows()
    ]

    t0 = time.perf_counter()
    results = {}
    with ProcessPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(eval_rhizaeon_single, t): t[0] for t in tasks}
        done = 0
        for fut in as_completed(futures):
            res = fut.result()
            results[res["row_idx"]] = res
            done += 1
            if done % 400 == 0 or done == total:
                elapsed = time.perf_counter() - t0
                print(f"[{done}/{total}] ({done/total*100:.1f}%) Elapsed: {elapsed:.1f}s")

    # Update DataFrame
    for idx, res in results.items():
        df.loc[idx, "rhizaeon_detected"] = res["rhizaeon_detected"]
        df.loc[idx, "rhizaeon_num_bps"] = res["rhizaeon_num_bps"]
        df.loc[idx, "rhizaeon_bps"] = res["rhizaeon_bps"]
        df.loc[idx, "rhizaeon_runtime_ms"] = res["rhizaeon_runtime_ms"]

    df.to_csv(csv_path, index=False)
    print(f"[*] Successfully updated: {csv_path}")

    # Recompute summary
    summary_rows = []
    scenarios = build_posada_scenarios()
    for sc in scenarios:
        sub = df[df["scenario_id"] == sc.scenario_id]
        if len(sub) == 0:
            continue
        summary_rows.append({
            "scenario_id": sc.scenario_id,
            "category": sc.category,
            "theta": sc.theta,
            "rho": sc.rho,
            "alpha": sc.alpha if sc.alpha is not None else "inf",
            "n_reps": len(sub),
            "mean_true_partitions": round(sub["true_num_partitions"].mean(), 2),
            "rhizaeon_detection_pct": round(sub["rhizaeon_detected"].mean() * 100.0, 1),
            "three_seq_detection_pct": round(sub["three_seq_detected"].mean() * 100.0, 1),
            "rhizaeon_mean_time_ms": round(sub["rhizaeon_runtime_ms"].mean(), 2),
            "three_seq_mean_time_ms": round(sub["three_seq_runtime_ms"].mean(), 2)
        })

    sum_df = pd.DataFrame(summary_rows)
    sum_csv = Path("benchmarks/posada_2001_benchmark/posada_summary.csv")
    sum_df.to_csv(sum_csv, index=False)
    print(f"[*] Successfully saved updated summary to: {sum_csv}")


if __name__ == "__main__":
    main()
