"""
benchmarks/posada_2001_benchmark/run_posada_benchmark.py
========================================================
High-throughput benchmark runner replicating Posada & Crandall (2001, PNAS):
Evaluates RhizAeon and 3SEQ across the 32 scenarios (Power & False Positive grids).
"""

import os
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
import argparse
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
from benchmarks.posada_2001_benchmark.generator import generate_posada_dataset, NT_CHARS
from benchmarks.posada_2001_benchmark.scenarios import build_posada_scenarios, PosadaScenario

THREE_SEQ_BIN = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq")
PTABLE_PATH = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq_ptable_250")
MASTER_SEED = 20011120 # Posada & Crandall PNAS publication date


def parse_3seq_output(csv_path: str) -> List[Dict[str, Any]]:
    """Parses 3SEQ output CSV file (run.3s.rec.csv)."""
    if not os.path.exists(csv_path):
        return []
    candidates = []
    try:
        with open(csv_path, "r") as f:
            lines = [l.strip() for l in f if l.strip()]
        if len(lines) <= 1:
            return []
        for line in lines[1:]:
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 7:
                p_val = float(parts[6]) if parts[6] != "" else 1.0
                candidates.append({
                    "P": parts[0],
                    "Q": parts[1],
                    "C": parts[2],
                    "p_val": p_val
                })
    except Exception:
        pass
    return candidates


def run_single_replicate(
    scenario: PosadaScenario,
    rep_idx: int,
    seed: int,
    run_3seq: bool = True
) -> Dict[str, Any]:
    """Runs a single replicate evaluation for both RhizAeon and 3SEQ."""
    # 1. Generate dataset
    t0_gen = time.perf_counter()
    mat, taxa, meta = generate_posada_dataset(
        n=scenario.n_taxa,
        L=scenario.length_nt,
        theta=scenario.theta,
        rho=scenario.rho,
        alpha=scenario.alpha,
        seed=seed
    )
    gen_time_ms = (time.perf_counter() - t0_gen) * 1000.0

    # 2. Run RhizAeon
    t0_rhiz = time.perf_counter()
    bps = []
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

    # Spatial accuracy calculation
    true_bps = meta.get("true_breakpoints", [])
    spatial_errs = []
    if meta["is_recombinant"] and len(true_bps) > 0 and len(rhiz_bps_coords) > 0:
        for t_bp in true_bps:
            min_e = min(abs(p_bp - t_bp) for p_bp in rhiz_bps_coords)
            spatial_errs.append(min_e)
    mean_spatial_err = float(np.mean(spatial_errs)) if spatial_errs else -1.0

    # 3. Run 3SEQ
    three_seq_detected = False
    three_seq_num_triplets = 0
    three_seq_time_ms = 0.0

    if run_3seq and os.path.exists(THREE_SEQ_BIN):
        t0_3s = time.perf_counter()
        try:
            with tempfile.TemporaryDirectory() as td:
                fa_path = os.path.join(td, "alignment.fa")
                with open(fa_path, "w") as f:
                    for i, name in enumerate(taxa):
                        seq_str = "".join(NT_CHARS[c] for c in mat[i])
                        f.write(f">{name}\n{seq_str}\n")
                        
                cmd = [
                    THREE_SEQ_BIN,
                    "-full", fa_path,
                    "-p", PTABLE_PATH,
                    "-id", f"rep_{rep_idx}",
                    "-q",
                    "-L40"
                ]
                subprocess.run(cmd, cwd=td, capture_output=True, text=True, check=False)
                three_seq_time_ms = (time.perf_counter() - t0_3s) * 1000.0
                
                csv_path = os.path.join(td, f"rep_{rep_idx}.3s.rec.csv")
                candidates = parse_3seq_output(csv_path)
                three_seq_detected = len(candidates) > 0
                three_seq_num_triplets = len(candidates)
        except Exception:
            three_seq_time_ms = (time.perf_counter() - t0_3s) * 1000.0

    return {
        "scenario_id": scenario.scenario_id,
        "category": scenario.category,
        "theta": scenario.theta,
        "rho": scenario.rho,
        "alpha": scenario.alpha if scenario.alpha is not None else "inf",
        "rep_idx": rep_idx,
        "seed": seed,
        "true_is_recombinant": meta["is_recombinant"],
        "true_num_partitions": meta["num_partitions"],
        "true_num_breakpoints": len(meta["true_breakpoints"]),
        "rhizaeon_detected": rhiz_detected,
        "rhizaeon_num_bps": rhiz_num_bps,
        "rhizaeon_bps": rhiz_bps_str,
        "rhizaeon_spatial_err_nt": mean_spatial_err,
        "rhizaeon_runtime_ms": rhiz_time_ms,
        "three_seq_detected": three_seq_detected,
        "three_seq_num_triplets": three_seq_num_triplets,
        "three_seq_runtime_ms": three_seq_time_ms,
        "gen_runtime_ms": gen_time_ms
    }


def main():
    parser = argparse.ArgumentParser(description="Replicate Posada & Crandall (2001) benchmark.")
    parser.add_argument("--workers", type=int, default=10, help="Number of concurrent workers.")
    parser.add_argument("--reps", type=int, default=100, help="Number of replicates per scenario.")
    parser.add_argument("--skip-3seq", action="store_true", help="Skip 3SEQ execution.")
    parser.add_argument("--limit-scenarios", type=int, default=None, help="Limit number of scenarios for quick testing.")
    parser.add_argument("--out-dir", type=str, default="benchmarks/posada_2001_benchmark", help="Output directory.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    scenarios = build_posada_scenarios()
    if args.limit_scenarios:
        scenarios = scenarios[:args.limit_scenarios]

    print(f"[*] Total scenarios: {len(scenarios)}, Replicates per scenario: {args.reps}")
    print(f"[*] Workers: {args.workers}, Running 3SEQ: {not args.skip_3seq}")

    tasks = []
    for s_idx, sc in enumerate(scenarios):
        for rep in range(args.reps):
            seed = MASTER_SEED + s_idx * 1000 + rep
            tasks.append((sc, rep, seed))

    print(f"[*] Total task instances: {len(tasks)}")
    t_start = time.perf_counter()

    results = []
    # Execute with concurrency <= 3
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(run_single_replicate, sc, rep, seed, not args.skip_3seq): (sc.scenario_id, rep)
            for sc, rep, seed in tasks
        }
        done_count = 0
        total_tasks = len(tasks)
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            done_count += 1
            if done_count % 100 == 0 or done_count == total_tasks:
                elapsed = time.perf_counter() - t_start
                rate = done_count / elapsed
                rem = (total_tasks - done_count) / rate if rate > 0 else 0
                print(f"[{done_count}/{total_tasks}] ({(done_count/total_tasks)*100:.1f}%) "
                      f"Elapsed: {elapsed:.1f}s, Est remaining: {rem:.1f}s")

    # Convert to DataFrame
    df = pd.DataFrame(results)
    raw_csv = out_dir / "posada_raw_results.csv"
    df.to_csv(raw_csv, index=False)
    print(f"[*] Saved raw results to: {raw_csv}")

    # Summary table
    summary_rows = []
    for sc in scenarios:
        sub = df[df["scenario_id"] == sc.scenario_id]
        if len(sub) == 0:
            continue
        rhiz_rate = sub["rhizaeon_detected"].mean() * 100.0
        three_rate = sub["three_seq_detected"].mean() * 100.0
        rhiz_time = sub["rhizaeon_runtime_ms"].mean()
        three_time = sub["three_seq_runtime_ms"].mean()
        mean_parts = sub["true_num_partitions"].mean()

        valid_errs = sub[sub["rhizaeon_spatial_err_nt"] >= 0]["rhizaeon_spatial_err_nt"]
        mean_spat_err = float(valid_errs.mean()) if len(valid_errs) > 0 else -1.0

        summary_rows.append({
            "scenario_id": sc.scenario_id,
            "category": sc.category,
            "theta": sc.theta,
            "rho": sc.rho,
            "alpha": sc.alpha if sc.alpha is not None else "inf",
            "n_reps": len(sub),
            "mean_true_partitions": round(mean_parts, 2),
            "rhizaeon_detection_pct": round(rhiz_rate, 1),
            "three_seq_detection_pct": round(three_rate, 1),
            "rhizaeon_mean_spatial_err_nt": round(mean_spat_err, 2),
            "rhizaeon_mean_time_ms": round(rhiz_time, 2),
            "three_seq_mean_time_ms": round(three_time, 2)
        })

    sum_df = pd.DataFrame(summary_rows)
    sum_csv = out_dir / "posada_summary.csv"
    sum_df.to_csv(sum_csv, index=False)
    print(f"[*] Saved summary to: {sum_csv}")


if __name__ == "__main__":
    main()
