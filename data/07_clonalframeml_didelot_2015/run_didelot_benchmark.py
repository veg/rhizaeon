"""
benchmarks/didelot_2015_benchmark/run_didelot_benchmark.py
==========================================================
Systematic replication of the ClonalFrame / ClonalFrameML / SimBac
bacterial microevolution and homologous gene conversion benchmark:
Didelot & Wilson (2015, PLoS Comput Biol 11:e1004041)
Didelot & Falush (2007, Genetics 175:1251-1266)

Evaluates:
1. Recombination Intensity: r/m in {0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0} (delta=1000, nu=0.02).
2. Import Tract Length: delta in {250, 500, 1000, 2500, 5000} (r/m=1.0, nu=0.02).
3. Import Divergence: nu in {0.005, 0.01, 0.03, 0.08} (r/m=1.0, delta=1000).
4. Taxonomic Cohort Scaling: N in {10, 20, 50, 100} on 25 kb bacterial genomes.
"""

import os
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
from pathlib import Path
from typing import Dict, Any, Tuple, List
from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
from benchmarks.didelot_2015_benchmark.generator import simulate_bacterial_gene_conversion


def evaluate_single_bacterial_genome(args: Tuple[str, str, int, int, int, float, float, float, float, int]) -> Dict[str, Any]:
    """Evaluates a single bacterial genome alignment under RhizAeon."""
    regime, cell_id, rep_id, n, L, r_over_m, delta, nu, theta, seed = args

    # 1. Simulate bacterial sequence evolution with gene conversion
    seqs, tracts, meta = simulate_bacterial_gene_conversion(
        n=n, L=L, r_over_m=r_over_m, delta=delta, nu=nu, theta=theta, seed=seed
    )
    true_bps = meta["true_breakpoints"]
    is_recombinant = meta["is_recombinant"]
    num_imports = meta["num_imports"]

    # 2. Modern RhizAeon
    char_map = {"A": 0, "C": 1, "G": 2, "T": 3}
    mat = np.array([[char_map.get(c, 0) for c in s] for s in seqs], dtype=np.int8)
    taxa_names = [f"Isolate_{i}" for i in range(n)]

    t0 = time.perf_counter()
    engine = PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=True)
    bps = run_recursive_partition_fda_screen(
        engine,
        taxa_names=taxa_names,
        min_z=1.8,
        min_pir=0.06,
        crossover_validation=True,
        crossover_p_threshold=0.005,
        min_informative_sites=3,
        polish_ml=True
    )
    rhiz_time_ms = (time.perf_counter() - t0) * 1000.0
    detected = len(bps) > 0
    num_detected_bps = len(bps)
    detected_coords = [b.polished_bp for b in bps]

    # Spatial accuracy evaluation
    if is_recombinant and len(detected_coords) > 0 and len(true_bps) > 0:
        errors = [min(abs(d - t) for t in true_bps) for d in detected_coords]
        min_error = min(errors)
        mean_error = float(np.mean(errors))
        within_50 = sum(e <= 50 for e in errors) / len(errors) * 100.0
        within_100 = sum(e <= 100 for e in errors) / len(errors) * 100.0
        within_250 = sum(e <= 250 for e in errors) / len(errors) * 100.0
    else:
        min_error = np.nan
        mean_error = np.nan
        within_50 = np.nan
        within_100 = np.nan
        within_250 = np.nan

    return {
        "regime": regime,
        "cell_id": cell_id,
        "rep_id": rep_id,
        "n_taxa": n,
        "genome_length": L,
        "r_over_m": r_over_m,
        "delta": delta,
        "nu": nu,
        "theta": theta,
        "is_recombinant": is_recombinant,
        "num_true_imports": num_imports,
        "num_true_bps": len(true_bps),
        "rhiz_detected": detected,
        "rhiz_num_bps": num_detected_bps,
        "rhiz_bps": ";".join(map(str, detected_coords)),
        "rhiz_min_spatial_error": min_error,
        "rhiz_mean_spatial_error": mean_error,
        "rhiz_acc_within_50nt": within_50,
        "rhiz_acc_within_100nt": within_100,
        "rhiz_acc_within_250nt": within_250,
        "rhiz_time_ms": rhiz_time_ms
    }


def main():
    bench_dir = Path(__file__).resolve().parent
    print("================================================================================")
    print("STARTING DIDELOT & WILSON (2015) BACTERIAL GENE CONVERSION BENCHMARK")
    print("================================================================================")

    tasks = []
    task_id = 1

    # Regime 1: Recombination Intensity r/m in {0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0}
    # (delta=1000, nu=0.02, theta=0.01, N=20, L=25000, 15 reps per cell = 105 datasets)
    for r_over_m in [0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0]:
        cell = f"intensity_rm_{r_over_m}"
        for rep in range(15):
            seed = 400000 + task_id
            tasks.append(("recombination_intensity", cell, rep, 20, 25000, r_over_m, 1000.0, 0.02, 0.01, seed))
            task_id += 1

    # Regime 2: Import Tract Length delta in {250, 500, 1000, 2500, 5000}
    # (r/m=1.0, nu=0.02, theta=0.01, N=20, L=25000, 15 reps per cell = 75 datasets)
    for delta in [250.0, 500.0, 1000.0, 2500.0, 5000.0]:
        cell = f"tract_length_d_{int(delta)}"
        for rep in range(15):
            seed = 500000 + task_id
            tasks.append(("tract_length", cell, rep, 20, 25000, 1.0, delta, 0.02, 0.01, seed))
            task_id += 1

    # Regime 3: Import Divergence nu in {0.005, 0.01, 0.03, 0.08}
    # (r/m=1.0, delta=1000, theta=0.01, N=20, L=25000, 15 reps per cell = 60 datasets)
    for nu in [0.005, 0.01, 0.03, 0.08]:
        cell = f"divergence_nu_{nu}"
        for rep in range(15):
            seed = 600000 + task_id
            tasks.append(("import_divergence", cell, rep, 20, 25000, 1.0, 1000.0, nu, 0.01, seed))
            task_id += 1

    # Regime 4: Taxonomic Cohort Scaling N in {10, 20, 50, 100}
    # (r/m=1.0, delta=1000, nu=0.02, theta=0.01, L=25000, 5 reps per cell = 20 datasets)
    for n in [10, 20, 50, 100]:
        cell = f"taxon_scaling_n_{n}"
        for rep in range(5):
            seed = 700000 + task_id
            tasks.append(("taxon_scaling", cell, rep, n, 25000, 1.0, 1000.0, 0.02, 0.01, seed))
            task_id += 1

    print(f"Total bacterial benchmark tasks generated: {len(tasks)}")
    print(f"  - Recombination intensity (r/m in [0..10]): 105 datasets")
    print(f"  - Tract length variation (delta in [250..5000]): 75 datasets")
    print(f"  - Import divergence (nu in [0.005..0.08]): 60 datasets")
    print(f"  - Taxon scaling (N in [10..100]): 20 datasets")

    # Parallel execution
    t_start = time.time()
    results = []
    max_workers = min(os.cpu_count() or 4, 8)
    print(f"Executing across {max_workers} worker processes...")

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(evaluate_single_bacterial_genome, t): t for t in tasks}
        completed_cnt = 0
        for f in as_completed(futures):
            res = f.result()
            results.append(res)
            completed_cnt += 1
            if completed_cnt % 25 == 0 or completed_cnt == len(tasks):
                elapsed = time.time() - t_start
                rate = completed_cnt / elapsed
                print(f"  Processed {completed_cnt}/{len(tasks)} ({completed_cnt/len(tasks)*100:.1f}%) in {elapsed:.1f}s ({rate:.2f} genomes/s)")

    # Save raw results
    raw_df = pd.DataFrame(results)
    raw_csv = bench_dir / "didelot_raw_results.csv"
    raw_df.to_csv(raw_csv, index=False)
    print(f"\nSaved raw results ({len(raw_df)} rows) to: {raw_csv}")

    # Compute cell summaries
    summary_rows = []
    for cell, g in raw_df.groupby("cell_id"):
        regime = g["regime"].iloc[0]
        n_taxa = g["n_taxa"].iloc[0]
        r_over_m = g["r_over_m"].iloc[0]
        delta = g["delta"].iloc[0]
        nu = g["nu"].iloc[0]
        n_reps = len(g)
        n_rec = g["is_recombinant"].sum()
        n_clonal = n_reps - n_rec

        tpr = (g[g["is_recombinant"]]["rhiz_detected"].sum() / n_rec * 100.0) if n_rec > 0 else np.nan
        fpr = (g[~g["is_recombinant"]]["rhiz_detected"].sum() / n_clonal * 100.0) if n_clonal > 0 else np.nan
        mean_min_mae = g["rhiz_min_spatial_error"].dropna().mean()
        mean_avg_mae = g["rhiz_mean_spatial_error"].dropna().mean()
        acc_50 = g["rhiz_acc_within_50nt"].dropna().mean()
        acc_100 = g["rhiz_acc_within_100nt"].dropna().mean()
        acc_250 = g["rhiz_acc_within_250nt"].dropna().mean()
        mean_time = g["rhiz_time_ms"].mean()
        mean_bps = g["rhiz_num_bps"].mean()

        summary_rows.append({
            "regime": regime,
            "cell_id": cell,
            "n_taxa": n_taxa,
            "r_over_m": r_over_m,
            "delta": delta,
            "nu": nu,
            "num_reps": n_reps,
            "num_recombinant": n_rec,
            "num_clonal": n_clonal,
            "rhiz_sensitivity": round(tpr, 2) if not np.isnan(tpr) else np.nan,
            "rhiz_fpr": round(fpr, 2) if not np.isnan(fpr) else np.nan,
            "rhiz_mean_bps": round(mean_bps, 2),
            "rhiz_min_spatial_error": round(mean_min_mae, 2) if not np.isnan(mean_min_mae) else np.nan,
            "rhiz_mean_spatial_error": round(mean_avg_mae, 2) if not np.isnan(mean_avg_mae) else np.nan,
            "rhiz_acc_within_50nt": round(acc_50, 2) if not np.isnan(acc_50) else np.nan,
            "rhiz_acc_within_100nt": round(acc_100, 2) if not np.isnan(acc_100) else np.nan,
            "rhiz_acc_within_250nt": round(acc_250, 2) if not np.isnan(acc_250) else np.nan,
            "rhiz_time_ms": round(mean_time, 2)
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_csv = bench_dir / "didelot_summary.csv"
    summary_df.to_csv(summary_csv, index=False)
    print(f"Saved summary ({len(summary_df)} cells) to: {summary_csv}")
    print(f"Total benchmark finished in {time.time() - t_start:.2f} seconds.")


if __name__ == "__main__":
    main()
