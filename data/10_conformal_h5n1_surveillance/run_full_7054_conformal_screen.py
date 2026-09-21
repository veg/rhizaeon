#!/usr/bin/env python3
"""
Full Panzootic Conformal Prediction Screen across All 7,054 Complete Genomes.

Applies distribution-free non-conformity scoring calibrated on N_calib = 400
to detect novel segment introductions, outgroup spillovers, and evolutionary outliers
across the entire global H5N1 panzootic database.
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

ALIGN_DIR = "/Users/sergei/Dropbox/Work/Collaborations/PittH5N1/2024-08-27/results/alignments"
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
LENS = [2277, 2271, 2148, 1701, 1494, 1407, 756, 690]

BENCH_DIR = REPO_ROOT / "benchmarks/flu_reassortment_benchmark"
RES_DIR = BENCH_DIR / "results"
FIG_DIR = REPO_ROOT / "figures"
RES_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 95)
print("  FULL-COHORT CONFORMAL PREDICTION SCREEN (7,054 COMPLETE 8-SEGMENT GENOMES)")
print("=" * 95)

# Step 1: Ingest all 8 segment alignments
print("\n[*] Ingesting 8 segment alignments from PittH5N1 archive...")
t0 = time.perf_counter()
seg_indices = {}
for s in SEGMENTS:
    fpath = os.path.join(ALIGN_DIR, f"{s}.fas")
    d = {}
    for r in SeqIO.parse(fpath, "fasta"):
        parts = r.id.split("|")
        if len(parts) >= 3:
            st = parts[2].strip()
            d[st] = str(r.seq).upper()
    seg_indices[s] = d

common_strains = sorted(list(set.intersection(*[set(d.keys()) for d in seg_indices.values()])))
N_total = len(common_strains)
print(f"[✓] Successfully ingested {N_total:,} complete 8-segment genomes in {time.perf_counter()-t0:.2f} s")

# Categorization helper
def categorize_strain(st):
    s = st.lower()
    if "texas/37" in s or "michigan/90" in s or "human" in s:
        return "Human"
    if "bear" in s:
        return "Bear"
    if "seal" in s or "sea_lion" in s:
        return "Marine_Mammal"
    if "cat" in s or "feline" in s:
        return "Cat"
    if "skunk" in s:
        return "Skunk"
    if "fox" in s:
        return "Fox"
    if "dairy_cow" in s or "cattle" in s or "bovine" in s:
        return "Bovine"
    if "eagle" in s or "falcon" in s or "hawk" in s or "owl" in s:
        return "Raptor"
    if "gull" in s or "tern" in s or "sandpiper" in s or "pelican" in s:
        return "Seabird_Shorebird"
    if "teal" in s or "wigeon" in s or "mallard" in s or "brant" in s or "duck" in s or "goose" in s or "swan" in s:
        return "Waterfowl"
    if "chicken" in s or "turkey" in s or "pheasant" in s or "quail" in s:
        return "Poultry"
    return "Other_Avian_Mammal"

strain_hosts = {st: categorize_strain(st) for st in common_strains}

# Step 2: Stratified Selection of Calibration Core (N_calib = 400)
print("\n[*] Assembling stratified calibration core (N_calib = 400)...")
np.random.seed(42)
host_buckets = {}
for st, h in strain_hosts.items():
    host_buckets.setdefault(h, []).append(st)

calib_strains = []
# Ensure all host clades are represented proportionally
for h, strains in host_buckets.items():
    n_sample = max(2, int(400 * (len(strains) / N_total)))
    n_sample = min(n_sample, len(strains))
    calib_strains.extend(np.random.choice(strains, size=n_sample, replace=False).tolist())

# Top up to exactly 400 if needed
remaining = [s for s in common_strains if s not in calib_strains]
if len(calib_strains) < 400:
    top_up = np.random.choice(remaining, size=400 - len(calib_strains), replace=False).tolist()
    calib_strains.extend(top_up)
else:
    calib_strains = calib_strains[:400]

print(f"[✓] Calibration core established with {len(calib_strains)} genomes across {len(host_buckets)} host groups.")

# Step 3: Discover Segment Reference Medoids on Calibration Core
print("\n[*] Stage 1 (Regime A): Discovering segment reference medoids on calibration core...")
t0 = time.perf_counter()
discovered_medoids = {}
calibration_non_conformity = {}

for s_idx, s_name in enumerate(SEGMENTS):
    # Encode calibration sequences as byte arrays
    calib_bytes = np.array([
        np.frombuffer(seg_indices[s_name][st].encode("ascii"), dtype=np.uint8)
        for st in calib_strains
    ])
    n_c, l_s = calib_bytes.shape
    
    # Compute pairwise distance matrix on calibration core
    # Vectorized Hamming p-distance
    D_calib = np.zeros((n_c, n_c), dtype=float)
    for i in range(n_c):
        D_calib[i] = np.mean(calib_bytes != calib_bytes[i], axis=1)
        
    # Unsupervised clustering: adaptive K (2 to 5)
    best_k = 2
    best_sil = -1.0
    best_labels = None
    for k in range(2, 6):
        clust = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average")
        lbls = clust.fit_predict(D_calib)
        sil = float(silhouette_score(D_calib, lbls, metric="precomputed"))
        if sil > best_sil:
            best_sil = sil
            best_k = k
            best_labels = lbls
            
    # Find exemplar medoids
    medoid_strains = []
    medoid_seq_arrays = []
    for c in range(best_k):
        members = np.where(best_labels == c)[0]
        sub_D = D_calib[np.ix_(members, members)]
        best_local = int(np.argmin(sub_D.sum(axis=1)))
        m_idx = members[best_local]
        medoid_strains.append(calib_strains[m_idx])
        medoid_seq_arrays.append(calib_bytes[m_idx])
        
    # Compute non-conformity scores on calibration set
    # Distance from each calibration sequence to nearest medoid
    dists_to_medoids = np.zeros((n_c, best_k), dtype=float)
    for c_id, m_arr in enumerate(medoid_seq_arrays):
        dists_to_medoids[:, c_id] = np.mean(calib_bytes != m_arr, axis=1)
        
    non_conf = np.sort(np.min(dists_to_medoids, axis=1))
    
    discovered_medoids[s_name] = {
        "k": best_k,
        "silhouette": best_sil,
        "strains": medoid_strains,
        "arrays": np.array(medoid_seq_arrays) # shape (K, L_s)
    }
    calibration_non_conformity[s_name] = non_conf
    
    conformal_threshold_95 = float(np.percentile(non_conf, 95.0))
    conformal_threshold_99 = float(np.percentile(non_conf, 99.0))
    print(f"    • {s_name}: K = {best_k} medoids (Sil: {best_sil:.3f}, Threshold α=0.05: {conformal_threshold_95:.4f}, α=0.01: {conformal_threshold_99:.4f})")

print(f"[✓] Reference medoid discovery complete in {time.perf_counter()-t0:.2f} s")

# Step 4: Stream Conformal Inference across All 7,054 Genomes
print(f"\n[*] Stage 2 (Regime B): Streaming Conformal Scoring across all {N_total:,} genomes...")
t0 = time.perf_counter()

# We evaluate segment by segment to bound memory usage
segment_p_values = np.zeros((N_total, len(SEGMENTS)), dtype=float)
segment_min_dists = np.zeros((N_total, len(SEGMENTS)), dtype=float)
segment_closest_medoid = np.zeros((N_total, len(SEGMENTS)), dtype=int)

for s_idx, s_name in enumerate(SEGMENTS):
    medoid_arr = discovered_medoids[s_name]["arrays"] # shape (K, L_s)
    calib_scores = calibration_non_conformity[s_name]
    n_calib = len(calib_scores)
    K = medoid_arr.shape[0]
    
    # Process queries in batches of 1,000 to maximize cache locality
    batch_size = 1000
    for b_start in range(0, N_total, batch_size):
        b_end = min(b_start + batch_size, N_total)
        batch_strains = common_strains[b_start:b_end]
        
        batch_bytes = np.array([
            np.frombuffer(seg_indices[s_name][st].encode("ascii"), dtype=np.uint8)
            for st in batch_strains
        ]) # shape (B, L_s)
        
        # Distances to K medoids: shape (B, K)
        dists_b = np.zeros((b_end - b_start, K), dtype=float)
        for k_idx in range(K):
            dists_b[:, k_idx] = np.mean(batch_bytes != medoid_arr[k_idx], axis=1)
            
        min_d = np.min(dists_b, axis=1)
        closest_m = np.argmin(dists_b, axis=1)
        
        # Conformal p-values: (1 + count(calib >= query)) / (N_calib + 1)
        # searchsorted gives count < query
        counts_less = np.searchsorted(calib_scores, min_d, side="left")
        counts_geq = n_calib - counts_less
        p_vals = (1.0 + counts_geq) / (n_calib + 1.0)
        
        segment_p_values[b_start:b_end, s_idx] = p_vals
        segment_min_dists[b_start:b_end, s_idx] = min_d
        segment_closest_medoid[b_start:b_end, s_idx] = closest_m

t_scoring = time.perf_counter() - t0
print(f"[✓] Conformal scoring completed across all 56,432 segments in {t_scoring:.2f} s ({t_scoring/N_total*1000:.3f} ms/genome)!")

# Step 5: Omnibus Hypothesis Aggregation
print("\n[*] Performing Simes and Bonferroni omnibus hypothesis aggregation...")
sorted_p = np.sort(segment_p_values, axis=1)
k_factors = 8.0 / np.arange(1, 9, dtype=float)
simes_matrix = sorted_p * k_factors[np.newaxis, :]
p_simes = np.clip(np.min(simes_matrix, axis=1), 0.0, 1.0)
p_bonf = np.clip(np.min(segment_p_values, axis=1) * 8.0, 0.0, 1.0)
p_min = np.min(segment_p_values, axis=1)

# Build comprehensive DataFrame
records = []
alpha_bonf = 0.05 / 8.0 # 0.00625

for i in range(N_total):
    st = common_strains[i]
    host = strain_hosts[st]
    
    outlier_segs = [
        SEGMENTS[s_idx] for s_idx in range(8)
        if segment_p_values[i, s_idx] <= alpha_bonf
    ]
    
    constellation = "-".join([f"M{segment_closest_medoid[i, s]+1}" for s in range(8)])
    
    rec = {
        "Strain": st,
        "Host": host,
        "Constellation": constellation,
        "Omnibus_Simes_P": float(p_simes[i]),
        "Omnibus_Bonf_P": float(p_bonf[i]),
        "Min_Segment_P": float(p_min[i]),
        "Is_Novel_FDR05": bool(p_simes[i] < 0.05),
        "Is_Novel_FDR01": bool(p_simes[i] < 0.01),
        "Novel_Segments": ",".join(outlier_segs) if outlier_segs else "None",
        "Num_Novel_Segments": len(outlier_segs)
    }
    for s_idx, s_name in enumerate(SEGMENTS):
        rec[f"{s_name}_P"] = float(segment_p_values[i, s_idx])
        rec[f"{s_name}_Dist"] = float(segment_min_dists[i, s_idx])
    records.append(rec)

df_all = pd.DataFrame(records)
catalog_path = RES_DIR / "full_7054_conformal_catalog.tsv"
df_all.to_csv(catalog_path, sep="\t", index=False)
print(f"[✓] Full 7,054 Conformal Catalog saved to {catalog_path}")

# Step 6: Detailed Statistical Audit of Findings
n_novel_05 = int(df_all["Is_Novel_FDR05"].sum())
n_novel_01 = int(df_all["Is_Novel_FDR01"].sum())

print("\n" + "=" * 90)
print("  FULL-COHORT CONFORMAL SURVEILLANCE FINDINGS")
print("=" * 90)
print(f"Total Isolates Screened: {N_total:,}")
print(f"Novel Genomic Outliers at FDR α = 0.05: {n_novel_05:,} ({n_novel_05/N_total*100:.2f}%)")
print(f"Novel Genomic Outliers at FDR α = 0.01: {n_novel_01:,} ({n_novel_01/N_total*100:.2f}%)")

print("\nNovelty Breakdown by Host Clade (FDR α = 0.05):")
host_novelty = df_all.groupby("Host")["Is_Novel_FDR05"].agg(["count", "sum", "mean"]).reset_index()
host_novelty.columns = ["Host", "Total", "Novel_Count", "Novel_Rate"]
host_novelty["Novel_Rate_%"] = (host_novelty["Novel_Rate"] * 100).round(2)
host_novelty = host_novelty.sort_values("Novel_Count", ascending=False)
print(host_novelty[["Host", "Total", "Novel_Count", "Novel_Rate_%"]].to_string(index=False))

print("\nNovel Segment Incursion Frequency across the Genome:")
seg_counts = {}
for s in SEGMENTS:
    # Segment is an outlier if its p-value <= alpha_bonf
    cnt = int((df_all[f"{s}_P"] <= alpha_bonf).sum())
    seg_counts[s] = cnt
    print(f"  • Segment {s:>3} ({LENS[SEGMENTS.index(s)]} nt): {cnt:,} isolates carry novel alleles ({cnt/N_total*100:.2f}%)")

print("\nTop 10 Extreme Outlier Genomes (Smallest Conformal Omnibus P-Values):")
top_outliers = df_all.sort_values(["Min_Segment_P", "Omnibus_Simes_P"]).head(10)
for idx, r in top_outliers.iterrows():
    print(f"  Strain: {r['Strain']}")
    print(f"    Host: {r['Host']} | Simes p = {r['Omnibus_Simes_P']:.4e} | Flagged Segments: [{r['Novel_Segments']}]")
    max_d_seg = max([(r[f'{s}_Dist'], s) for s in SEGMENTS])
    print(f"    Peak Divergence: {max_d_seg[1]} (distance = {max_d_seg[0]:.4f})")

# Step 7: Generate Publication-Grade Figure
print("\n[*] Generating Full-Cohort Conformal Surveillance Figure...")
fig, axes = plt.subplots(1, 3, figsize=(19, 5.8))

# Panel A: Omnibus Simes P-value Distribution across Host Clades
ax = axes[0]
major_hosts = ["Poultry", "Waterfowl", "Bovine", "Raptor", "Marine_Mammal", "Cat", "Fox", "Skunk"]
host_p_data = [df_all[df_all["Host"] == h]["Omnibus_Simes_P"].values for h in major_hosts]
bp = ax.boxplot(host_p_data, labels=major_hosts, patch_artist=True, medianprops=dict(color="black", linewidth=1.5))
colors = ["#457b9d", "#2a9d8f", "#e76f51", "#f4a261", "#9b5de5", "#e63946", "#8338ec", "#3a86ff"]
for patch, col in zip(bp["boxes"], colors):
    patch.set_facecolor(col)
    patch.set_alpha(0.7)
ax.axhline(0.05, color="red", linestyle="--", linewidth=1.5, label=r"Novelty Cutoff ($\alpha=0.05$)")
ax.set_xticklabels(major_hosts, rotation=35, ha="right")
ax.set_ylabel("Omnibus Simes Conformal p-value", fontsize=11)
ax.set_title("(A) Genome Conformal Validity across Hosts", fontsize=12)
ax.legend(loc="upper right")
ax.grid(axis="y", linestyle="--", alpha=0.4)

# Panel B: Segment Novelty Burden
ax = axes[1]
x_pos = np.arange(len(SEGMENTS))
counts_list = [seg_counts[s] for s in SEGMENTS]
rates_list = [c / N_total * 100 for c in counts_list]
bars = ax.bar(x_pos, rates_list, color="#2b5c8f", alpha=0.85, edgecolor="black")
ax.set_xticks(x_pos)
ax.set_xticklabels(SEGMENTS)
ax.set_ylabel("Isolates with Novel Allele (%)", fontsize=11)
ax.set_title("(B) Reassortant Novelty Burden per Segment", fontsize=12)
ax.grid(axis="y", linestyle="--", alpha=0.4)
for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2.0, yval + 0.1, f"{yval:.1f}%", ha='center', va='bottom', fontsize=9)

# Panel C: Multi-Segment Reassortant Burden (Number of Novel Segments per Isolate)
ax = axes[2]
burden_counts = df_all["Num_Novel_Segments"].value_counts().sort_index()
burden_df = pd.DataFrame({"Segments": range(9)})
burden_df["Count"] = burden_df["Segments"].map(burden_counts).fillna(0)
burden_df["Pct"] = burden_df["Count"] / N_total * 100

bars_c = ax.bar(burden_df["Segments"], burden_df["Pct"], color="#e76f51", alpha=0.85, edgecolor="black")
ax.set_yscale("log")
ax.set_xlabel("Number of Statistically Novel Segments", fontsize=11)
ax.set_ylabel("Isolates (%) [Log Scale]", fontsize=11)
ax.set_title("(C) Panzootic Reassortant Incursion Depth", fontsize=12)
ax.grid(axis="y", linestyle="--", alpha=0.4)
for bar, pct in zip(bars_c, burden_df["Pct"]):
    if pct > 0:
        ax.text(bar.get_x() + bar.get_width()/2.0, pct * 1.2, f"{pct:.1f}%", ha='center', va='bottom', fontsize=8)

plt.tight_layout()
fig_png = FIG_DIR / "fig_full_7054_conformal_novelty.png"
fig_pdf = FIG_DIR / "fig_full_7054_conformal_novelty.pdf"
plt.savefig(fig_png, dpi=300)
plt.savefig(fig_pdf)
plt.close()
print(f"[✓] Generated publication figures: {fig_png.name} and {fig_pdf.name}")

print("\n[✓] Full-Cohort Conformal Prediction Screen successfully completed!")
