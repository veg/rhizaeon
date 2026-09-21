#!/usr/bin/env python3
"""
benchmarks/grand_1000_benchmark/run_5_paradigms_benchmark.py
============================================================
Benchmark Harness comparing the 5 Algorithmic Detection Paradigms of RhizAeon
across the Grand-Scale 1,000-Alignment Ground-Truth Codon Benchmark suite:
  1. Classical MDS + Orthogonal Procrustes + L-PIR (Baseline RhizAeonDetector)
  2. Two-Tier Cascaded (Fast Scalar Tier 1 -> Static 384D Embedding Tier 2)
  3. Recursive Partitioning FDA (RP-FDA / 1D Total Variation)
  4. Continuous Multiresolution Wavelet (Haar-Procrustes WTMM + Ridge Refinement)
  5. Adaptive 3-Tier Cascade (FDA Triage -> Wavelet Zoom -> Procrustes Refinement)
  [6. 3SEQ as external benchmark reference]

Evaluates 1,000 alignments (20 scenarios x 50 replicates):
  - True Power (Sensitivity on 800 recombinant alignments)
  - Empirical False Positive Rate (FPR on 200 null confounder alignments)
  - Localization Error (MAE in codons, exact, <=1 cod, <=3 cod)
  - Recombinant Taxon Identification Accuracy
  - Execution Latency (ms)
"""

import os
import sys
import json
import time
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from rhizaeon.tensor import parse_fasta
from rhizaeon.embed import build_prefix_engine
from rhizaeon.segmentation import RhizAeonDetector
from rhizaeon.fda import run_recursive_partition_fda_screen
from rhizaeon.wavelet import run_wavelet_recombination_screen
from rhizaeon.adaptive import run_adaptive_hybrid_screen
from rhizaeon.pir import refine_breakpoint_codon

BENCH_DIR = REPO_ROOT / "benchmarks" / "grand_1000_benchmark"
DATA_DIR = BENCH_DIR / "data"
RES_DIR = BENCH_DIR / "results"
RES_DIR.mkdir(parents=True, exist_ok=True)

RAW_CSV = RES_DIR / "grand_1000_5_paradigms_raw.csv"
SUMMARY_CSV = RES_DIR / "grand_1000_5_paradigms_summary.csv"
DIM_SUMMARY_CSV = RES_DIR / "grand_1000_5_paradigms_dimension_summary.csv"
LOG_FILE = RES_DIR / "benchmark_5_paradigms.log"

EXISTING_1000_RAW = RES_DIR / "grand_1000_raw_results.csv"


def log(msg: str):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line, flush=True)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")


def compute_metrics(detected_bps: List[int], detected_taxa: List[str], true_bps: List[int]) -> Dict[str, Any]:
    """Computes rigorous quantitative detection accuracy metrics."""
    has_det = len(detected_bps) > 0
    is_null = len(true_bps) == 0

    if is_null:
        return {
            "has_detection": has_det,
            "false_alarm": has_det,
            "min_err": np.nan,
            "mean_err": np.nan,
            "within_0": not has_det,
            "within_1": not has_det,
            "within_3": not has_det,
            "within_15": not has_det,
            "within_30": not has_det,
            "recombinant_correct": not has_det
        }

    if not has_det:
        return {
            "has_detection": False,
            "false_alarm": False,
            "min_err": np.nan,
            "mean_err": np.nan,
            "within_0": False,
            "within_1": False,
            "within_3": False,
            "within_15": False,
            "within_30": False,
            "recombinant_correct": False
        }

    errs = []
    for tb in true_bps:
        dists = [abs(tb - db) for db in detected_bps]
        errs.append(min(dists))

    min_err = float(np.min(errs))
    mean_err = float(np.mean(errs))
    rec_correct = any("R_" in dt or "mosaic" in dt.lower() for dt in detected_taxa)

    return {
        "has_detection": True,
        "false_alarm": False,
        "min_err": min_err,
        "mean_err": mean_err,
        "within_0": min_err == 0,
        "within_1": min_err <= 1,
        "within_3": min_err <= 3,
        "within_15": min_err <= 15,
        "within_30": min_err <= 30,
        "recombinant_correct": rec_correct
    }


def eval_mds(fa_path: str, taxa: List[str], detector: RhizAeonDetector) -> Tuple[List[int], List[str], float]:
    t0 = time.perf_counter()
    eng = build_prefix_engine(fa_path, engine="scalar", codon=True)
    events = detector.detect_recombination(eng, taxa)
    runtime = time.perf_counter() - t0
    detected_bps = sorted(list(set([e["breakpoint"] for e in events])))
    detected_taxa = sorted(list(set([e["recombinant"] for e in events])))
    return detected_bps, detected_taxa, runtime


def eval_twotier(fa_path: str, taxa: List[str], detector: RhizAeonDetector) -> Tuple[List[int], List[str], float]:
    t0 = time.perf_counter()
    eng = build_prefix_engine(fa_path, engine="two-tier", codon=True)
    events = detector.detect_recombination(eng, taxa)
    runtime = time.perf_counter() - t0
    detected_bps = sorted(list(set([e["breakpoint"] for e in events])))
    detected_taxa = sorted(list(set([e["recombinant"] for e in events])))
    return detected_bps, detected_taxa, runtime


def eval_fda(fa_path: str, taxa: List[str]) -> Tuple[List[int], List[str], float]:
    t0 = time.perf_counter()
    eng = build_prefix_engine(fa_path, engine="scalar", codon=True)
    bps = run_recursive_partition_fda_screen(eng, taxa, min_z=2.0, min_pir=0.20, min_flank=25)
    runtime = time.perf_counter() - t0
    detected_bps = sorted(list(set([b.breakpoint_nt for b in bps])))
    detected_taxa = sorted(list(set([b.recombinant_taxon for b in bps])))
    return detected_bps, detected_taxa, runtime


def eval_wavelet(fa_path: str, taxa: List[str]) -> Tuple[List[int], List[str], float]:
    t0 = time.perf_counter()
    eng = build_prefix_engine(fa_path, engine="scalar", codon=True)
    scalo, ridges = run_wavelet_recombination_screen(
        eng, taxa, min_scale=32, max_scale=256, step=15, min_prominence=0.5, min_z_threshold=2.5, pir_threshold=0.20
    )
    # Require composite score >= 1.5
    strong_ridges = [r for r in ridges if r.l_pir >= 0.20 and (r.max_z_score * r.l_pir) >= 1.5]
    detected_bps = []
    detected_taxa = []
    for r in strong_ridges[:3]:
        tax_idx = taxa.index(r.taxon_name)
        p1_idx = taxa.index(r.parent_1) if r.parent_1 in taxa else 0
        p2_idx = taxa.index(r.parent_2) if r.parent_2 in taxa else 1
        ref_bp, ref_pir = refine_breakpoint_codon(eng, r.singularity_nt, tax_idx, p1_idx, p2_idx, search_radius=20, flank_len=25)
        if ref_pir >= 0.20:
            detected_bps.append(ref_bp)
            detected_taxa.append(r.taxon_name)
    runtime = time.perf_counter() - t0
    return sorted(list(set(detected_bps))), sorted(list(set(detected_taxa))), runtime


def eval_adaptive(fa_path: str, taxa: List[str]) -> Tuple[List[int], List[str], float]:
    t0 = time.perf_counter()
    eng = build_prefix_engine(fa_path, engine="scalar", codon=True)
    bps = run_adaptive_hybrid_screen(eng, taxa, fda_bin_size=80, fda_step=15, min_kinetic_z=2.0, min_pir=0.20)
    runtime = time.perf_counter() - t0
    detected_bps = sorted(list(set([b.breakpoint_nt for b in bps])))
    detected_taxa = sorted(list(set([b.recombinant_taxon for b in bps])))
    return detected_bps, detected_taxa, runtime


def main():
    manifest_path = DATA_DIR / "benchmark_manifest.json"
    if not manifest_path.exists():
        log(f"Error: manifest not found at {manifest_path}")
        sys.exit(1)

    with open(manifest_path, "r") as f:
        manifest = json.load(f)

    log("=" * 95)
    log(f"  BENCHMARKING 5 RHIZAEON PARADIGMS ON {len(manifest):,} GROUND-TRUTH ALIGNMENTS")
    log("=" * 95)

    # Checkpoint loading
    completed_keys = set()
    raw_rows = []
    if RAW_CSV.exists():
        try:
            df_exist = pd.read_csv(RAW_CSV)
            raw_rows = df_exist.to_dict("records")
            for r in raw_rows:
                completed_keys.add((str(r["scenario"]), int(r["replicate"]), str(r["method"])))
            log(f"[✓] Checkpoint loaded: {len(completed_keys)} evaluations already in {RAW_CSV.name}")
        except Exception as e:
            log(f"Warning loading checkpoint: {e}")

    # Import existing Classical MDS and 3SEQ from grand_1000_raw_results.csv if available
    if EXISTING_1000_RAW.exists():
        try:
            df_orig = pd.read_csv(EXISTING_1000_RAW)
            for _, r in df_orig.iterrows():
                m_name = str(r["method"])
                if m_name == "rhizaeon_scalar":
                    target_method = "1_classical_mds"
                elif m_name == "three_seq":
                    target_method = "6_three_seq"
                else:
                    continue
                key = (str(r["scenario"]), int(r["replicate"]), target_method)
                if key not in completed_keys:
                    rec = r.to_dict()
                    rec["method"] = target_method
                    raw_rows.append(rec)
                    completed_keys.add(key)
            log(f"[✓] Integrated existing baseline rows. Total checkpointed: {len(completed_keys)}")
        except Exception as e:
            log(f"Warning importing existing rows: {e}")

    detector = RhizAeonDetector(
        window_units=25,
        min_tract_units=30,
        ghost_z_threshold=2.75,
        pir_threshold=0.20,
        k_dims=4
    )

    t_start = time.perf_counter()
    save_interval = 20

    for d_idx, ds in enumerate(manifest):
        scen = ds["scenario"]
        rep = ds["replicate"]
        dim = ds["dimension"]
        fa_path = ds["alignment_path"]
        true_bps = ds.get("true_bps_codon", [])
        is_null = (len(true_bps) == 0)

        taxa, _ = parse_fasta(fa_path)

        methods_to_run = [
            ("1_classical_mds", lambda: eval_mds(fa_path, taxa, detector)),
            ("2_two_tier", lambda: eval_twotier(fa_path, taxa, detector)),
            ("3_rp_fda", lambda: eval_fda(fa_path, taxa)),
            ("4_wavelet_wtmm", lambda: eval_wavelet(fa_path, taxa)),
            ("5_adaptive_3tier", lambda: eval_adaptive(fa_path, taxa))
        ]

        for method_name, fn in methods_to_run:
            key = (scen, rep, method_name)
            if key in completed_keys:
                continue

            try:
                det_bps, det_taxa, runtime = fn()
            except Exception as e:
                log(f"[!] Error in {method_name} on {scen} rep {rep}: {e}")
                det_bps, det_taxa, runtime = [], [], 0.0

            metrics = compute_metrics(det_bps, det_taxa, true_bps)
            row = {
                "dataset_id": ds["dataset_id"],
                "scenario": scen,
                "dimension": dim,
                "replicate": rep,
                "method": method_name,
                "is_null": is_null,
                "num_true_bps": len(true_bps),
                "true_bps": str(true_bps),
                "detected_bps": str(det_bps),
                "detected_taxa": str(det_taxa),
                "runtime_sec": round(runtime, 5),
                **metrics
            }
            raw_rows.append(row)
            completed_keys.add(key)

        if (d_idx + 1) % save_interval == 0 or (d_idx + 1) == len(manifest):
            pd.DataFrame(raw_rows).to_csv(RAW_CSV, index=False)
            elapsed = time.perf_counter() - t_start
            pct = 100.0 * (d_idx + 1) / len(manifest)
            log(f"[{pct:5.1f}%] Processed {d_idx+1}/{len(manifest)} datasets | Total evaluated rows: {len(raw_rows)} | Elapsed: {elapsed:.1f}s")

    # Final save
    df_all = pd.DataFrame(raw_rows)
    df_all.to_csv(RAW_CSV, index=False)
    log(f"[✓] Raw results successfully saved to {RAW_CSV}")

    # =========================================================================
    # COMPUTE SUMMARY AGGREGATIONS
    # =========================================================================
    summary_rows = []
    methods_order = ["1_classical_mds", "2_two_tier", "3_rp_fda", "4_wavelet_wtmm", "5_adaptive_3tier", "6_three_seq"]
    display_names = {
        "1_classical_mds": "Classical MDS (Scalar Tier 1)",
        "2_two_tier": "Two-Tier Cascaded (Hybrid)",
        "3_rp_fda": "RP-FDA (1D Total Variation)",
        "4_wavelet_wtmm": "Continuous Wavelets (WTMM)",
        "5_adaptive_3tier": "Adaptive 3-Tier Cascade",
        "6_three_seq": "3SEQ (Exhaustive Triplet)"
    }

    for m in methods_order:
        df_m = df_all[df_all["method"] == m]
        if len(df_m) == 0:
            continue

        df_null = df_m[df_m["is_null"] == True]
        df_rec = df_m[df_m["is_null"] == False]

        fpr = (df_null["false_alarm"].mean() * 100.0) if len(df_null) > 0 else np.nan
        power = (df_rec["has_detection"].mean() * 100.0) if len(df_rec) > 0 else np.nan

        mae = df_rec[df_rec["has_detection"]]["min_err"].mean()
        exact = (df_rec["within_0"].mean() * 100.0) if len(df_rec) > 0 else np.nan
        w1 = (df_rec["within_1"].mean() * 100.0) if len(df_rec) > 0 else np.nan
        w3 = (df_rec["within_3"].mean() * 100.0) if len(df_rec) > 0 else np.nan
        tax_acc = (df_rec["recombinant_correct"].mean() * 100.0) if len(df_rec) > 0 else np.nan
        latency_ms = df_m["runtime_sec"].mean() * 1000.0

        summary_rows.append({
            "Method Key": m,
            "Architecture / Paradigm": display_names.get(m, m),
            "Power (%)": round(power, 1),
            "FPR (%)": round(fpr, 1),
            "MAE (codons)": round(mae, 2) if not np.isnan(mae) else np.nan,
            "Exact (%)": round(exact, 1),
            "<=1 Codon (%)": round(w1, 1),
            "<=3 Codons (%)": round(w3, 1),
            "Taxon Acc (%)": round(tax_acc, 1),
            "Latency (ms)": round(latency_ms, 1)
        })

    df_summary = pd.DataFrame(summary_rows)
    df_summary.to_csv(SUMMARY_CSV, index=False)
    log(f"\n{'='*95}\n  OVERALL SUMMARY (1,000 ALIGNMENTS)\n{'='*95}")
    log("\n" + df_summary.to_string(index=False))

    # =========================================================================
    # GENERATE LATEX TABLE
    # =========================================================================
    tex_path = REPO_ROOT / "paper" / "tables" / "tab_5_paradigms_1000_benchmark.tex"
    tex_path.parent.mkdir(parents=True, exist_ok=True)
    with open(tex_path, "w") as f:
        f.write("% Systematic evaluation of 5 RhizAeon paradigms across 1,000 ground-truth codon alignments\n")
        f.write("\\begin{table*}[t!]\n\\centering\n\\footnotesize\n")
        f.write("{\\setlength{\\tabcolsep}{3.5pt}\n")
        f.write("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}l cccccccc@{}}\n")
        f.write("\\toprule\n")
        f.write("Paradigm & Power (\\%) & FPR (\\%) & MAE (cod) & Exact & $\\le$1 Cod & $\\le$3 Cod & Taxon Acc & Latency \\\\\n")
        f.write("\\midrule\n")
        for _, r in df_summary.iterrows():
            f.write(f"{r['Architecture / Paradigm']} & {r['Power (%)']}\\% & {r['FPR (%)']}\\% & {r['MAE (codons)']:.2f} & {r['Exact (%)']}\\% & {r['<=1 Codon (%)']}\\% & {r['<=3 Codons (%)']}\\% & {r['Taxon Acc (%)']}\\% & {r['Latency (ms)']:.1f} ms \\\\\n")
        f.write("\\bottomrule\n")
        f.write("\\end{tabular*}\n}\n")
        f.write("\\caption{Systematic evaluation of five \\RhizAeon algorithmic detection paradigms across the Grand-Scale 1,000-Alignment Recombination Benchmark suite. ")
        f.write("Evaluations span twenty evolutionary scenarios (four null confounder tiers and sixteen recombinant conditions) with fifty replicates each (1,000 alignments total; 16 taxa, 600 codons). ")
        f.write("Power denotes detection sensitivity on recombinant datasets. FPR denotes false alarm rate on negative controls under spatial rate variation (Gamma $\\alpha=0.3$) and selection heterogeneity ($\\omega=0.85$). ")
        f.write("MAE measures mean absolute error to the closest true breakpoint in codons. Latency reports mean CPU execution time per complete alignment.}\n")
        f.write("\\label{tab:5_paradigms_1000_benchmark}\n")
        f.write("\\end{table*}\n")
    log(f"[✓] LaTeX table saved to {tex_path}")


if __name__ == "__main__":
    main()
