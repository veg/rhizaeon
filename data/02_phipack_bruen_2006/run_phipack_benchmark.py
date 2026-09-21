"""
simulations/02_phipack_bruen_2006/run_phipack_benchmark.py
==========================================================
High-throughput benchmark runner replicating Bruen et al. (2006, Genetics):
Evaluates RhizAeon, PHI (Normal), MaxChi^2, NSS, and 3SEQ across 2,300 alignments
spanning coalescent recombination power, neutral constant nulls, exponential growth,
and Gamma rate heterogeneity.
"""

import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"

import sys
import time
import argparse
import tempfile
import re
import subprocess
from pathlib import Path
from typing import Dict, Any, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from rhizaeon.tensor import PrefixDistanceEngine
from rhizaeon.fda import run_recursive_partition_fda_screen
from generator import generate_phipack_dataset, NT_CHARS

PHI_BIN = "/Users/sergei/miniconda3/bin/Phi"
THREE_SEQ_BIN = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq")
PTABLE_PATH = str(PROJECT_ROOT / "benchmarks" / "tools" / "3seq" / "cmake-build" / "3seq_ptable_250")
MASTER_SEED = 20060801  # Bruen et al. Genetics August 2006 publication date


def build_phipack_scenarios(reps_per_condition: int = 50) -> List[Dict[str, Any]]:
    """Builds the complete experimental grid of 2,300 alignment conditions."""
    scenarios = []

    # 1. Power Grid: n in {10, 50}, theta in {5, 10, 20}, rho in {1, 4, 16, 64}
    for n in [10, 50]:
        for theta in [5.0, 10.0, 20.0]:
            for rho in [1.0, 4.0, 16.0, 64.0]:
                for r in range(reps_per_condition):
                    scenarios.append({
                        "scenario_id": f"power_n{n}_th{int(theta)}_rho{int(rho)}",
                        "category": "power",
                        "n_taxa": n,
                        "length_nt": 1000,
                        "theta": theta,
                        "rho": rho,
                        "beta": 0.0,
                        "alpha": None,
                        "rep_idx": r
                    })

    # 2. Constant Neutral Null: n in {10, 50}, theta in {5, 10, 20}, rho = 0
    for n in [10, 50]:
        for theta in [5.0, 10.0, 20.0]:
            for r in range(reps_per_condition):
                scenarios.append({
                    "scenario_id": f"null_const_n{n}_th{int(theta)}",
                    "category": "null_constant",
                    "n_taxa": n,
                    "length_nt": 1000,
                    "theta": theta,
                    "rho": 0.0,
                    "beta": 0.0,
                    "alpha": None,
                    "rep_idx": r
                })

    # 3. Exponential Growth Null: n in {10, 50}, theta in {10, 20}, rho = 0, beta in {5.0, 20.0}
    for n in [10, 50]:
        for theta in [10.0, 20.0]:
            for beta in [5.0, 20.0]:
                for r in range(reps_per_condition):
                    scenarios.append({
                        "scenario_id": f"null_growth_n{n}_th{int(theta)}_beta{int(beta)}",
                        "category": "null_growth",
                        "n_taxa": n,
                        "length_nt": 1000,
                        "theta": theta,
                        "rho": 0.0,
                        "beta": beta,
                        "alpha": None,
                        "rep_idx": r
                    })

    # 4. Gamma Rate Heterogeneity Null: n in {10, 50}, theta in {10, 20}, rho = 0, alpha in {0.5, 2.0}
    for n in [10, 50]:
        for theta in [10.0, 20.0]:
            for alpha in [0.5, 2.0]:
                for r in range(reps_per_condition):
                    scenarios.append({
                        "scenario_id": f"null_gamma_n{n}_th{int(theta)}_alpha{alpha}",
                        "category": "null_gamma",
                        "n_taxa": n,
                        "length_nt": 1000,
                        "theta": theta,
                        "rho": 0.0,
                        "beta": 0.0,
                        "alpha": alpha,
                        "rep_idx": r
                    })

    return scenarios


def evaluate_single_replicate(item: Tuple[int, Dict[str, Any]]) -> Dict[str, Any]:
    """Runs a single simulation replicate across all 5 benchmarked methods."""
    global_idx, sc = item
    seed = (MASTER_SEED + global_idx * 7919) % (2**31 - 1)

    t0 = time.perf_counter()

    # 1. Generate alignment via msprime
    msa, taxa, meta = generate_phipack_dataset(
        n=sc["n_taxa"],
        L=sc["length_nt"],
        theta=sc["theta"],
        rho=sc["rho"],
        beta=sc["beta"],
        alpha=sc["alpha"],
        seed=seed
    )

    # 2. RhizAeon
    t_rhiz0 = time.perf_counter()
    rhiz_detected = False
    rhiz_num_bps = 0
    rhiz_bps_str = ""
    try:
        engine = PrefixDistanceEngine(msa, codon_aligned=False, compute_transitions=True)
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
        rhiz_detected = len(bps) > 0
        rhiz_num_bps = len(bps)
        coords = [b.polished_bp if b.polished_bp is not None else b.breakpoint_nt for b in bps]
        rhiz_bps_str = ";".join(str(c) for c in coords)
    except Exception:
        pass
    rhiz_time_ms = (time.perf_counter() - t_rhiz0) * 1000.0

    # Write temporary FASTA for PhiPack and 3SEQ inside dedicated temp dir
    phi_p = 1.0
    maxchi_p = 1.0
    nss_p = 1.0
    phi_time_ms = 0.0

    three_p = 1.0
    three_detected = False
    three_time_ms = 0.0

    with tempfile.TemporaryDirectory() as td:
        fa_path = os.path.join(td, "alignment.fasta")
        with open(fa_path, "w") as f:
            for i, name in enumerate(taxa):
                f.write(f">{name}\n" + "".join(NT_CHARS[c] for c in msa[i]) + "\n")

        # 3. PhiPack (Phi, MaxChi, NSS)
        t_phi0 = time.perf_counter()
        try:
            res_phi = subprocess.run(
                [PHI_BIN, "-f", fa_path, "-t", "D", "-o"],
                cwd=td,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120
            )
            out_phi = res_phi.stdout
            m_phi = re.search(r"PHI \(Normal\):\s+([0-9.eE+-]+)", out_phi)
            m_maxchi = re.search(r"Max Chi\^2:\s+([0-9.eE+-]+)", out_phi)
            m_nss = re.search(r"NSS:\s+([0-9.eE+-]+)", out_phi)

            if m_phi:
                phi_p = float(m_phi.group(1))
            if m_maxchi:
                maxchi_p = float(m_maxchi.group(1))
            if m_nss:
                nss_p = float(m_nss.group(1))
        except Exception:
            pass
        phi_time_ms = (time.perf_counter() - t_phi0) * 1000.0

        # 4. 3SEQ
        t_three0 = time.perf_counter()
        try:
            run_id = f"r{global_idx}"
            subprocess.run(
                [THREE_SEQ_BIN, "-full", fa_path, "-p", PTABLE_PATH, "-id", run_id, "-q", "-L40"],
                cwd=td,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=120
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
                            three_p = min(three_p, p_val)
                            if p_val < 0.05:
                                three_detected = True
        except Exception:
            pass
        three_time_ms = (time.perf_counter() - t_three0) * 1000.0

    total_time_ms = (time.perf_counter() - t0) * 1000.0

    return {
        "global_idx": global_idx,
        "scenario_id": sc["scenario_id"],
        "category": sc["category"],
        "rep_idx": sc["rep_idx"],
        "seed": seed,
        "n_taxa": sc["n_taxa"],
        "length_nt": sc["length_nt"],
        "theta": sc["theta"],
        "rho": sc["rho"],
        "beta": sc["beta"],
        "alpha": sc["alpha"] if sc["alpha"] is not None else "inf",
        "is_recombinant": meta["is_recombinant"],
        "num_true_trees": meta["num_trees"],
        "num_true_bps": meta["num_bps"],
        "segregating_sites": meta["segregating_sites"],
        "informative_sites": meta["informative_sites"],
        # RhizAeon
        "rhiz_detected": rhiz_detected,
        "rhiz_num_bps": rhiz_num_bps,
        "rhiz_bps": rhiz_bps_str,
        "rhiz_time_ms": round(rhiz_time_ms, 2),
        # PHI
        "phi_p_val": phi_p,
        "phi_detected": (phi_p < 0.05),
        # MaxChi
        "maxchi_p_val": maxchi_p,
        "maxchi_detected": (maxchi_p < 0.05),
        # NSS
        "nss_p_val": nss_p,
        "nss_detected": (nss_p < 0.05),
        "phi_suite_time_ms": round(phi_time_ms, 2),
        # 3SEQ
        "three_p_val": three_p,
        "three_detected": three_detected,
        "three_time_ms": round(three_time_ms, 2),
        # Overall
        "total_time_ms": round(total_time_ms, 2)
    }


def main():
    parser = argparse.ArgumentParser(description="Bruen et al. 2006 PhiPack Benchmark Runner")
    parser.add_argument("--workers", type=int, default=8, help="Number of parallel worker processes")
    parser.add_argument("--reps", type=int, default=50, help="Replicates per parameter scenario")
    parser.add_argument("--out_dir", type=str, default=str(Path(__file__).resolve().parent), help="Output directory")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_csv = out_dir / "phipack_raw_results.csv"
    summary_csv = out_dir / "phipack_summary.csv"

    scenarios = build_phipack_scenarios(reps_per_condition=args.reps)
    print(f"Total alignment scenarios queued: {len(scenarios)} using {args.workers} workers")

    items = list(enumerate(scenarios))
    results = []
    t_start = time.perf_counter()

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(evaluate_single_replicate, it): it[0] for it in items}
        done_count = 0
        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)
            done_count += 1
            if done_count % 100 == 0 or done_count == len(scenarios):
                elapsed = time.perf_counter() - t_start
                rate = done_count / elapsed
                rem = (len(scenarios) - done_count) / max(1e-5, rate)
                print(f"[{done_count}/{len(scenarios)}] Completed in {elapsed:.1f}s ({rate:.1f} align/s, ~{rem:.1f}s remaining)")

    # Sort results
    results.sort(key=lambda x: x["global_idx"])
    df = pd.DataFrame(results)
    df.to_csv(raw_csv, index=False)
    print(f"Saved raw results to {raw_csv}")

    # Aggregated Summary by Category and Scenario
    summary_rows = []
    for (cat, n, th, rho, beta, alpha), grp in df.groupby(["category", "n_taxa", "theta", "rho", "beta", "alpha"]):
        n_reps = len(grp)
        summary_rows.append({
            "category": cat,
            "n_taxa": n,
            "theta": th,
            "rho": rho,
            "beta": beta,
            "alpha": alpha,
            "n_reps": n_reps,
            "mean_informative_sites": round(grp["informative_sites"].mean(), 1),
            "rhiz_rate": round(grp["rhiz_detected"].mean() * 100.0, 1),
            "phi_rate": round(grp["phi_detected"].mean() * 100.0, 1),
            "maxchi_rate": round(grp["maxchi_detected"].mean() * 100.0, 1),
            "nss_rate": round(grp["nss_detected"].mean() * 100.0, 1),
            "three_rate": round(grp["three_detected"].mean() * 100.0, 1),
            "mean_rhiz_ms": round(grp["rhiz_time_ms"].mean(), 2),
            "mean_phi_ms": round(grp["phi_suite_time_ms"].mean(), 2),
            "mean_three_ms": round(grp["three_time_ms"].mean(), 2)
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.sort_values(by=["category", "n_taxa", "theta", "rho", "beta", "alpha"], inplace=True)
    df_summary.to_csv(summary_csv, index=False)
    print(f"Saved aggregated summary to {summary_csv}")
    print("\n--- Summary Preview ---")
    print(df_summary.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
