#!/usr/bin/env python3
"""
benchmarks/grand_1000_benchmark/generate_grand_1000_data.py
===========================================================
Generates 1,000 ground-truth codon alignments across 20 evolutionary scenarios
(50 stochastic replicates each) using mechanistic pyvolve codon models.
"""

import os
import sys
import json
import time
import tempfile
import argparse
from pathlib import Path
from typing import Dict, List, Any
import pyvolve
from Bio import SeqIO

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from benchmarks.grand_1000_benchmark.scenarios import GRAND_SCENARIOS, scale_newick

DATA_DIR = REPO_ROOT / "benchmarks" / "grand_1000_benchmark" / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def simulate_single_dataset(scen_name: str, rep_idx: int, cfg: Dict[str, Any], out_dir: Path) -> Dict[str, Any]:
    rep_dir = out_dir / scen_name / f"rep_{rep_idx}"
    rep_dir.mkdir(parents=True, exist_ok=True)
    fa_path = rep_dir / "alignment.fasta"
    gt_path = rep_dir / "ground_truth.json"

    # If already generated, return existing metadata
    if fa_path.exists() and gt_path.exists():
        with open(gt_path, "r") as f:
            return json.load(f)

    tree_scale = cfg["tree_scale"]
    srv_alpha = cfg["srv_alpha"]
    omega = cfg["omega"]
    prune_taxa = set(cfg.get("prune_taxa", []))
    partitions_cfg = cfg["partitions"]

    taxa_seqs: Dict[str, str] = {}
    partition_meta = []
    curr_codon = 0

    with tempfile.TemporaryDirectory() as tmp_dir:
        for p_idx, (p_len, base_tree_str) in enumerate(partitions_cfg):
            scaled_tree_str = scale_newick(base_tree_str, tree_scale)
            py_tree = pyvolve.read_tree(tree=scaled_tree_str)

            mod_params = {"omega": omega, "kappa": 2.5}
            if srv_alpha is not None:
                mod_params["alpha"] = float(srv_alpha)
                mod_params["k"] = 4

            model = pyvolve.Model("codon", mod_params)
            part = pyvolve.Partition(models=model, size=p_len)
            evolver = pyvolve.Evolver(partitions=part, tree=py_tree)

            tmp_fa = os.path.join(tmp_dir, f"part_{p_idx}.fasta")
            evolver(seqfile=tmp_fa, count=1, quiet=True)

            for r in SeqIO.parse(tmp_fa, "fasta"):
                t_name = r.id
                if t_name in prune_taxa:
                    continue
                taxa_seqs[t_name] = taxa_seqs.get(t_name, "") + str(r.seq)

            partition_meta.append({
                "partition_idx": p_idx,
                "start_codon": curr_codon + 1,
                "end_codon": curr_codon + p_len,
                "length_codons": p_len
            })
            curr_codon += p_len

    # Write alignment
    with open(fa_path, "w") as f:
        for t_name in sorted(taxa_seqs.keys()):
            f.write(f">{t_name}\n{taxa_seqs[t_name]}\n")

    total_codons = curr_codon
    total_nt = total_codons * 3

    meta = {
        "dataset_id": f"{scen_name}_rep_{rep_idx}",
        "scenario": scen_name,
        "dimension": cfg["dimension"],
        "description": cfg["description"],
        "replicate": rep_idx,
        "num_taxa": len(taxa_seqs),
        "total_codons": total_codons,
        "total_nt": total_nt,
        "num_true_bps": len(cfg["true_bps"]),
        "true_bps_codon": cfg["true_bps"],
        "recombinant_taxa": cfg["recombinant_taxa"],
        "tree_scale": tree_scale,
        "srv_alpha": srv_alpha,
        "omega": omega,
        "alignment_path": str(fa_path),
        "partitions": partition_meta
    }

    with open(gt_path, "w") as f:
        json.dump(meta, f, indent=2)

    return meta


def main():
    parser = argparse.ArgumentParser(description="Generate 1,000 Ground-Truth Datasets")
    parser.add_argument("--reps_per_scenario", type=int, default=50, help="Replicates per scenario (default: 50 -> 1000 total)")
    args = parser.parse_args()

    n_scen = len(GRAND_SCENARIOS)
    total_target = n_scen * args.reps_per_scenario
    print("=" * 95)
    print(f"  GENERATING GRAND-SCALE BENCHMARK DATASET: {total_target:,} ALIGNMENTS")
    print(f"  {n_scen} Scenarios x {args.reps_per_scenario} Replicates")
    print("=" * 95)

    manifest_records = []
    t0 = time.perf_counter()

    for s_idx, (scen_name, cfg) in enumerate(GRAND_SCENARIOS.items()):
        t_scen = time.perf_counter()
        print(f"[{s_idx+1:02d}/{n_scen}] Generating '{scen_name}' ({args.reps_per_scenario} reps)...", end="", flush=True)
        for rep in range(args.reps_per_scenario):
            rec = simulate_single_dataset(scen_name, rep, cfg, DATA_DIR)
            manifest_records.append(rec)
        print(f" done in {time.perf_counter()-t_scen:.2f} s")

    elapsed = time.perf_counter() - t0
    print(f"\n[✓] Generated {len(manifest_records):,} datasets in {elapsed:.2f} s ({elapsed/len(manifest_records)*1000:.1f} ms/dataset)")

    # Save manifest
    manifest_json = DATA_DIR / "benchmark_manifest.json"
    manifest_csv = DATA_DIR / "benchmark_manifest.csv"
    with open(manifest_json, "w") as f:
        json.dump(manifest_records, f, indent=2)

    import pandas as pd
    pd.DataFrame(manifest_records).to_csv(manifest_csv, index=False)
    print(f"[✓] Saved manifest to {manifest_json} and {manifest_csv}")


if __name__ == "__main__":
    main()
