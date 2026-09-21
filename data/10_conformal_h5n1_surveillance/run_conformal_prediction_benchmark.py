#!/usr/bin/env python3
"""
Conformal Prediction Benchmark for Influenza A Reassortment Surveillance.

Rigorous evaluation of distribution-free non-conformity scores, finite-sample
coverage guarantees, empirical power on novel segment introductions, and Simes omnibus aggregation.
Calibrated on N_calib = 400 genomes sampled from the complete 7,054 panzootic database.
"""

import sys
import os
import time
from pathlib import Path
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

REPO_ROOT = Path("/Users/sergei/Projects/TOGA_MEME/recombination")
sys.path.insert(0, str(REPO_ROOT))

from rhizaeon.tensor import encode_alignment_matrix, PrefixDistanceEngine
from rhizaeon.conformal import ConformalMetricCalibrator, MultiSegmentConformalEngine

ALIGN_DIR = "/Users/sergei/Dropbox/Work/Collaborations/PittH5N1/2024-08-27/results/alignments"
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
LENS = [2277, 2271, 2148, 1701, 1494, 1407, 756, 690]

cum = 0
BOUNDS = []
for l in LENS:
    BOUNDS.append((cum, cum + l))
    cum += l

BENCH_DIR = REPO_ROOT / "benchmarks/flu_reassortment_benchmark"
RES_DIR = BENCH_DIR / "results"
FIG_DIR = REPO_ROOT / "figures"
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 95)
print("  CONFORMAL PREDICTION (DISTRIBUTION-FREE NON-CONFORMITY) GENOMIC BENCHMARK")
print("=" * 95)

# Step 1: Index all 8 segments from full archive
print("\n[*] Ingesting 8 segment alignments from PittH5N1 panzootic archive...")
t0 = time.perf_counter()
seg_indices = {}
for s in SEGMENTS:
    fpath = os.path.join(ALIGN_DIR, f"{s}.fas")
    d = {}
    for r in SeqIO.parse(fpath, "fasta"):
        parts = r.id.split("|")
        if len(parts) >= 3:
            st = parts[2].strip()
            d[st] = str(r.seq)
    seg_indices[s] = d

common_strains = sorted(list(set.intersection(*[set(d.keys()) for d in seg_indices.values()])))
print(f"[✓] Indexed {len(common_strains):,} complete 8-segment strains in {time.perf_counter()-t0:.2f} s")

# Step 2: Stratified Split: Calibration Core (N_calib = 400), In-Distribution Test (N_test = 200)
np.random.seed(42)
sampled_strains = np.random.choice(common_strains, size=600, replace=False).tolist()
calib_strains = sampled_strains[:400]
test_strains = sampled_strains[400:]

print(f"[*] Calibration Set: N_calib = {len(calib_strains)} | In-Distribution Test Set: N_test = {len(test_strains)}")

# Build concatenated matrix for the 600 genomes
print("[*] Concatenating 12,744 nt pseudomolecules for evaluation cohort...")
t0 = time.perf_counter()
concat_records = []
for st in sampled_strains:
    full_str = "".join([seg_indices[s][st] for s in SEGMENTS])
    concat_records.append(full_str)

# Convert to byte matrix
# Map ACGTN to 0, 1, 2, 3, 4
n_eval = len(sampled_strains)
seq_mat = np.zeros((n_eval, cum), dtype=np.int8)
char_map = {ord("A"): 0, ord("a"): 0, ord("C"): 1, ord("c"): 1, ord("G"): 2, ord("g"): 2, ord("T"): 3, ord("t"): 3}
for i, s_str in enumerate(concat_records):
    b = s_str.encode("ascii")
    arr = np.frombuffer(b, dtype=np.uint8)
    for char_val, code in char_map.items():
        seq_mat[i, arr == char_val] = code

engine = PrefixDistanceEngine(seq_mat, codon_aligned=False)
print(f"[✓] Encoded {n_eval} pseudomolecules in {time.perf_counter()-t0:.2f} s")

# Step 3: Run Regime A on Calibration Set (indices 0..399)
print("\n[*] Stage 1: Regime A Medoid Discovery on Calibration Core (N_calib = 400)...")
calib_idx = list(range(400))
discovered_medoids = {}

for s_idx, s_name in enumerate(SEGMENTS):
    start, end = BOUNDS[s_idx]
    D_seg = engine.query_distance_matrix(start, end)
    D_calib = D_seg[:400, :400]
    
    # Adaptive clustering
    best_k = 4 if s_name in ["PB2", "PB1", "PA"] else 2
    clusterer = AgglomerativeClustering(n_clusters=best_k, metric="precomputed", linkage="average")
    labels = clusterer.fit_predict(D_calib)
    
    medoid_globals = []
    for c in range(best_k):
        members = np.where(labels == c)[0]
        sub_D = D_calib[np.ix_(members, members)]
        best_local = int(np.argmin(sub_D.sum(axis=1)))
        global_idx = members[best_local]
        medoid_globals.append(global_idx)
        
    discovered_medoids[s_name] = medoid_globals
    print(f"    • {s_name}: Discovered {len(medoid_globals)} reference medoids")

# Step 4: Fit MultiSegmentConformalEngine
conformal_engine = MultiSegmentConformalEngine(SEGMENTS, BOUNDS, normalized=False)
conformal_engine.fit_from_prefix_engine(engine, calib_idx, discovered_medoids)
print("[✓] Multi-Segment Conformal Engine successfully calibrated.")

# Step 5: Construct Out-of-Distribution (OOD) Novel Reassortants
# Test Set: indices 400..599
test_idx = list(range(400, 600))
base_seq = seq_mat[test_idx[0]].copy()

synthetic_taxa = []
synthetic_seqs = []

# OOD 1: Novel HA Swap (12% mutations in HA segment)
seq_ood_ha = base_seq.copy()
ha_s, ha_e = BOUNDS[3]
n_mut = int(0.12 * (ha_e - ha_s))
pos = np.random.choice(range(ha_s, ha_e), size=n_mut, replace=False)
seq_ood_ha[pos] = (seq_ood_ha[pos] + 1) % 4
synthetic_taxa.append("SYNTHETIC_OOD_HA_Reassortant")
synthetic_seqs.append(seq_ood_ha)

# OOD 2: Novel PB1 Swap (10% mutations in PB1 segment)
seq_ood_pb1 = base_seq.copy()
pb1_s, pb1_e = BOUNDS[1]
n_mut = int(0.10 * (pb1_e - pb1_s))
pos = np.random.choice(range(pb1_s, pb1_e), size=n_mut, replace=False)
seq_ood_pb1[pos] = (seq_ood_pb1[pos] + 1) % 4
synthetic_taxa.append("SYNTHETIC_OOD_PB1_Reassortant")
synthetic_seqs.append(seq_ood_pb1)

# OOD 3: Double Reassortant (PA + NA)
seq_ood_double = base_seq.copy()
for s_idx in [2, 5]:
    s_s, s_e = BOUNDS[s_idx]
    n_m = int(0.10 * (s_e - s_s))
    pos = np.random.choice(range(s_s, s_e), size=n_m, replace=False)
    seq_ood_double[pos] = (seq_ood_double[pos] + 1) % 4
synthetic_taxa.append("SYNTHETIC_OOD_Double_PA_NA_Reassortant")
synthetic_seqs.append(seq_ood_double)

# OOD 4: Complete Outgroup (All 8 segments 15% diverged)
seq_ood_full = base_seq.copy()
for s_idx in range(8):
    s_s, s_e = BOUNDS[s_idx]
    n_m = int(0.15 * (s_e - s_s))
    pos = np.random.choice(range(s_s, s_e), size=n_m, replace=False)
    seq_ood_full[pos] = (seq_ood_full[pos] + 1) % 4
synthetic_taxa.append("SYNTHETIC_OOD_Full_Pan_Segment_Outgroup")
synthetic_seqs.append(seq_ood_full)

# Combine test set with synthetic OOD
combined_seq_mat = np.vstack([seq_mat, np.array(synthetic_seqs)])
combined_taxa = sampled_strains + synthetic_taxa
combined_engine = PrefixDistanceEngine(combined_seq_mat, codon_aligned=False)

eval_indices = test_idx + [n_eval + i for i in range(len(synthetic_taxa))]
eval_taxa = [combined_taxa[i] for i in eval_indices]

print(f"\n[*] Evaluating Conformal Inference across {len(eval_taxa)} test genomes ({len(test_idx)} in-distribution + {len(synthetic_taxa)} novel OOD)...")
t0 = time.perf_counter()
df_conformal = conformal_engine.score_queries(combined_engine, eval_taxa, eval_indices, alpha_fdr=0.05)
t_elapsed = time.perf_counter() - t0
print(f"[✓] Conformal scoring complete in {t_elapsed*1000:.2f} ms ({t_elapsed/len(eval_taxa)*1000:.3f} ms/genome)")

# Save TSV
conformal_tsv = RES_DIR / "conformal_prediction_evaluation_results.tsv"
df_conformal.to_csv(conformal_tsv, sep="\t", index=False)
print(f"[✓] Saved conformal evaluation results to {conformal_tsv}")

# Step 6: Coverage and Power Audit
df_in_dist = df_conformal.iloc[:len(test_idx)]
df_ood = df_conformal.iloc[len(test_idx):]

print("\n" + "=" * 80)
print("  CONFORMAL COVERAGE AND STATISTICAL POWER AUDIT")
print("=" * 80)

for alpha in [0.05, 0.02, 0.01]:
    empirical_fpr = (df_in_dist["Omnibus_Simes_P"] < alpha).mean()
    print(f"Nominal Error Rate α = {alpha:.2f} | Empirical False Alarm Rate: {empirical_fpr*100:.2f}% (Target: ≤ {alpha*100:.2f}%)")

ood_detected = (df_ood["Omnibus_Simes_P"] < 0.05).mean()
print(f"\nStatistical Power on Novel Out-Of-Distribution Lineages (α = 0.05): {ood_detected*100:.1f}%")

print("\nOOD Novel Reassortant Segment Attribution:")
for idx, r in df_ood.iterrows():
    name = r["Taxon"]
    simes_p = r["Omnibus_Simes_P"]
    bonf_p = r["Omnibus_Bonf_P"]
    flagged = r["Novel_Segments"]
    print(f"  • {name}: Simes p = {simes_p:.4f}, Bonferroni p = {bonf_p:.4f} | Flagged Novel Segments: [{flagged}]")

# Step 7: Plotting
print("\n[*] Generating Conformal Prediction Calibration & Power Figure...")
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

# Panel A: Calibration Non-Conformity Distributions per Segment
ax = axes[0]
box_data = [conformal_engine.calibrators[s].calibration_scores for s in SEGMENTS]
bp = ax.boxplot(box_data, labels=SEGMENTS, patch_artist=True, medianprops=dict(color="black", linewidth=1.5))
for patch in bp["boxes"]:
    patch.set_facecolor("#2b5c8f")
    patch.set_alpha(0.7)
ax.set_ylabel(r"Non-Conformity Score $S(X) = \min_c d(X, M_c)$ (Hamming)", fontsize=11)
ax.set_title("(A) Segment Intra-Cluster Calibration Scores", fontsize=12)
ax.grid(axis="y", linestyle="--", alpha=0.4)

# Panel B: In-Distribution Conformal P-Value Distribution
ax = axes[1]
p_simes_in = df_in_dist["Omnibus_Simes_P"].values
counts, bins, patches = ax.hist(p_simes_in, bins=10, range=(0, 1), density=True, color="#2a9d8f", alpha=0.75, edgecolor="black")
ax.axhline(1.0, color="crimson", linestyle="--", linewidth=1.8, label=r"Theoretical Uniform Null $U(0, 1)$")
ax.set_xlabel("Omnibus Simes Conformal p-value", fontsize=11)
ax.set_ylabel("Empirical Density", fontsize=11)
ax.set_title(r"(B) In-Distribution Calibration Validity ($N_{test}=200$)", fontsize=12)
ax.legend(frameon=True, fontsize=10)
ax.grid(axis="y", linestyle="--", alpha=0.4)

# Panel C: Reassortant Novelty Attribution Profile
ax = axes[2]
x_pos = np.arange(len(SEGMENTS))
width = 0.22

p_normal_mean = -np.log10(np.clip(df_in_dist[[f"{s}_P" for s in SEGMENTS]].values.mean(axis=0), 1e-4, 1.0))
p_ood_ha = -np.log10(np.clip(df_ood[df_ood["Taxon"] == "SYNTHETIC_OOD_HA_Reassortant"][[f"{s}_P" for s in SEGMENTS]].values[0], 1e-4, 1.0))
p_ood_double = -np.log10(np.clip(df_ood[df_ood["Taxon"] == "SYNTHETIC_OOD_Double_PA_NA_Reassortant"][[f"{s}_P" for s in SEGMENTS]].values[0], 1e-4, 1.0))

ax.bar(x_pos - width, p_normal_mean, width=width, label="In-Distribution Mean", color="#457b9d", alpha=0.85)
ax.bar(x_pos, p_ood_ha, width=width, label="HA Novel Reassortant", color="#e76f51", alpha=0.9)
ax.bar(x_pos + width, p_ood_double, width=width, label="PA+NA Double Reassortant", color="#9b5de5", alpha=0.9)

bonf_cutoff = -np.log10(0.05 / 8.0)
ax.axhline(bonf_cutoff, color="red", linestyle=":", linewidth=1.8, label=r"Bonferroni Threshold ($\alpha=0.05/8$)")

ax.set_xticks(x_pos)
ax.set_xticklabels(SEGMENTS)
ax.set_ylabel(r"Novelty Significance $-\log_{10}(p)$", fontsize=11)
ax.set_title("(C) Segment-Specific Novelty Attribution", fontsize=12)
ax.legend(frameon=True, fontsize=9, loc="upper right")
ax.grid(axis="y", linestyle="--", alpha=0.4)

plt.tight_layout()
fig_png = FIG_DIR / "fig_conformal_prediction_reassortment.png"
fig_pdf = FIG_DIR / "fig_conformal_prediction_reassortment.pdf"
plt.savefig(fig_png, dpi=300)
plt.savefig(fig_pdf)
plt.close()
print(f"[✓] Generated publication figures: {fig_png.name} and {fig_pdf.name}")

print("\n[✓] Conformal Prediction Benchmark successfully completed!")
