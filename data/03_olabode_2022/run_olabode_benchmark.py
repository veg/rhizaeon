"""
simulations/03_olabode_2022/run_olabode_benchmark.py
=====================================================
High-throughput benchmark runner executing the calibrated standard RhizAeon workflow
and 3SEQ across all 360 continuous-time HIV-1 datasets from Olabode et al. (2022, PNAS):
- Exp 1: 16-taxon reference benchmark (COMET design; 30 datasets)
- Exp 2: 37-taxon continuous time-scaled simulations with Gamma SRV (300 datasets)
- Exp 3: Taxon scaling cohorts (N=50, 100, 200 sequences; 30 datasets)
"""

import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import json
import time
import argparse
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from rhizaeon.tensor import PrefixDistanceEngine, encode_alignment_matrix
from rhizaeon.fda import run_recursive_partition_fda_screen
from eval_olabode_metrics import (
    compute_partition_misclassification_error,
    compute_breakpoint_spatial_metrics
)

THREE_SEQ_BIN = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq")
PTABLE_PATH = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq_ptable_250")


def evaluate_single_dataset(item: Tuple[int, Dict[str, Any]]) -> Dict[str, Any]:
    """Evaluates a single simulated HIV-1 alignment with RhizAeon and 3SEQ."""
    idx, entry = item
    fasta_path = entry["fasta_path"]
    true_bps_nt = entry.get("true_breakpoints_nt", [])
    num_taxa = entry["num_taxa"]
    length_nt = entry["length_nt"]

    # 1. Run RhizAeon Standard Calibrated Workflow
    t_rhiz0 = time.perf_counter()
    rhiz_bps_nt = []
    try:
        mat, taxa, L = encode_alignment_matrix(fasta_path)
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
        rhiz_bps_nt = [b.polished_bp if b.polished_bp is not None else b.breakpoint_nt for b in bps]
        rhiz_bps_nt = sorted(list(set(rhiz_bps_nt)))
    except Exception:
        pass
    rhiz_time_ms = (time.perf_counter() - t_rhiz0) * 1000.0

    rhiz_error_pct = compute_partition_misclassification_error(true_bps_nt, rhiz_bps_nt, length_nt)
    rhiz_spatial = compute_breakpoint_spatial_metrics(true_bps_nt, rhiz_bps_nt)

    # 2. Run 3SEQ
    t_three0 = time.perf_counter()
    three_bps_nt = []
    three_detected = False
    three_min_p = 1.0
    with tempfile.TemporaryDirectory() as td:
        run_id = f"ol_{idx}"
        try:
            subprocess.run(
                [THREE_SEQ_BIN, "-full", fasta_path, "-p", PTABLE_PATH, "-id", run_id, "-q", "-L40"],
                cwd=td,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=180
            )
            csv_path = os.path.join(td, f"{run_id}.3s.rec.csv")
            if os.path.exists(csv_path):
                with open(csv_path) as f:
                    lines = [l.strip() for l in f if l.strip()]
                if len(lines) > 1:
                    for l in lines[1:]:
                        parts = [p.strip() for p in l.split(",")]
                        if len(parts) >= 7 and parts[6] != "":
                            p_val = float(parts[6])
                            three_min_p = min(three_min_p, p_val)
                            if p_val < 0.05:
                                three_detected = True
                                # Add breakpoints (left and right)
                                if len(parts) >= 6:
                                    try:
                                        b1, b2 = int(parts[3]), int(parts[4])
                                        three_bps_nt.extend([b1, b2])
                                    except Exception:
                                        pass
        except Exception:
            pass
    three_time_ms = (time.perf_counter() - t_three0) * 1000.0

    three_bps_nt = sorted(list(set(three_bps_nt)))
    three_error_pct = compute_partition_misclassification_error(true_bps_nt, three_bps_nt, length_nt)
    three_spatial = compute_breakpoint_spatial_metrics(true_bps_nt, three_bps_nt)

    return {
        "global_idx": idx,
        "experiment": entry["experiment"],
        "scenario": entry.get("scenario", f"scaling_N{num_taxa}"),
        "replicate": entry["replicate"],
        "num_taxa": num_taxa,
        "length_nt": length_nt,
        "num_true_bps": len(true_bps_nt),
        "true_bps": ";".join(str(b) for b in true_bps_nt),
        # RhizAeon
        "rhiz_num_bps": len(rhiz_bps_nt),
        "rhiz_bps": ";".join(str(b) for b in rhiz_bps_nt),
        "rhiz_detected": len(rhiz_bps_nt) > 0,
        "rhiz_error_pct": round(rhiz_error_pct, 3),
        "rhiz_mean_spatial_error_nt": round(rhiz_spatial.get("mean_error_nt", np.nan), 1) if not np.isnan(rhiz_spatial.get("mean_error_nt", np.nan)) else None,
        "rhiz_median_spatial_error_nt": round(rhiz_spatial.get("median_error_nt", np.nan), 1) if not np.isnan(rhiz_spatial.get("median_error_nt", np.nan)) else None,
        "rhiz_time_ms": round(rhiz_time_ms, 2),
        # 3SEQ
        "three_num_bps": len(three_bps_nt),
        "three_bps": ";".join(str(b) for b in three_bps_nt),
        "three_detected": three_detected,
        "three_min_p": three_min_p,
        "three_error_pct": round(three_error_pct, 3),
        "three_mean_spatial_error_nt": round(three_spatial.get("mean_error_nt", np.nan), 1) if not np.isnan(three_spatial.get("mean_error_nt", np.nan)) else None,
        "three_time_ms": round(three_time_ms, 2)
    }


def main():
    parser = argparse.ArgumentParser(description="Olabode et al. (2022) HIV-1 Simulation Benchmark Runner")
    parser.add_argument("--workers", type=int, default=8, help="Number of worker processes")
    parser.add_argument("--manifest", type=str, default=str(Path(__file__).resolve().parent / "datasets_manifest.json"), help="Datasets manifest")
    parser.add_argument("--out_dir", type=str, default=str(Path(__file__).resolve().parent), help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = out_dir / "olabode_raw_results.csv"
    summary_csv = out_dir / "olabode_summary.csv"
    scaling_csv = out_dir / "olabode_scaling.csv"

    with open(args.manifest) as f:
        manifest_data = json.load(f)

    print(f"Loaded {len(manifest_data)} datasets from manifest. Evaluating across {args.workers} workers...")
    items = list(enumerate(manifest_data))
    results = []
    t_start = time.perf_counter()

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(evaluate_single_dataset, it): it[0] for it in items}
        done = 0
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            done += 1
            if done % 30 == 0 or done == len(manifest_data):
                elapsed = time.perf_counter() - t_start
                rate = done / elapsed
                rem = (len(manifest_data) - done) / max(1e-5, rate)
                print(f"[{done}/{len(manifest_data)}] Completed in {elapsed:.1f}s ({rate:.1f} align/s, ~{rem:.1f}s remaining)")

    results.sort(key=lambda x: x["global_idx"])
    df = pd.DataFrame(results)
    df.to_csv(raw_csv, index=False)
    print(f"Saved raw results to {raw_csv}")

    # Summary by experiment and scenario
    summary_rows = []
    for (exp, sc, n_bps, n_taxa), grp in df.groupby(["experiment", "scenario", "num_true_bps", "num_taxa"]):
        summary_rows.append({
            "experiment": exp,
            "scenario": sc,
            "num_true_bps": n_bps,
            "num_taxa": n_taxa,
            "n_reps": len(grp),
            "rhiz_detection_rate": round(grp["rhiz_detected"].mean() * 100.0, 1),
            "rhiz_median_error_pct": round(grp["rhiz_error_pct"].median(), 2),
            "rhiz_q25_error_pct": round(grp["rhiz_error_pct"].quantile(0.25), 2),
            "rhiz_q75_error_pct": round(grp["rhiz_error_pct"].quantile(0.75), 2),
            "rhiz_mean_time_ms": round(grp["rhiz_time_ms"].mean(), 1),
            "three_detection_rate": round(grp["three_detected"].mean() * 100.0, 1),
            "three_median_error_pct": round(grp["three_error_pct"].median(), 2),
            "three_mean_time_ms": round(grp["three_time_ms"].mean(), 1)
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.sort_values(by=["experiment", "num_true_bps", "num_taxa"], inplace=True)
    df_summary.to_csv(summary_csv, index=False)
    print(f"Saved summary to {summary_csv}")

    # Scaling table for Exp 3
    df_exp3 = df[df["experiment"] == "exp3_scaling"]
    scaling_rows = []
    for n_taxa, grp in df_exp3.groupby("num_taxa"):
        scaling_rows.append({
            "num_taxa": n_taxa,
            "n_reps": len(grp),
            "rhiz_time_ms": round(grp["rhiz_time_ms"].mean(), 1),
            "rhiz_time_sec": round(grp["rhiz_time_ms"].mean() / 1000.0, 3),
            "three_time_ms": round(grp["three_time_ms"].mean(), 1),
            "three_time_sec": round(grp["three_time_ms"].mean() / 1000.0, 3),
        })
    df_scale = pd.DataFrame(scaling_rows)
    if not df_scale.empty:
        df_scale.sort_values(by="num_taxa", inplace=True)
    df_scale.to_csv(scaling_csv, index=False)
    print(f"Saved scaling summary to {scaling_csv}")

    print("\n--- Summary Preview ---")
    print(df_summary.to_string(index=False))


if __name__ == "__main__":
    main()
