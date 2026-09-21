#!/usr/bin/env python3
"""
Test Dynamic Medoid Recruitment across Full 7,054 Influenza A Genomes.

Compares static reference typing against adaptive dynamic medoid recruitment:
- Evaluates alert rate reduction (combating alarm fatigue).
- Traces cluster foundation by true biological outgroups (e.g. Peruvian sea lions and ducks).
- Quantifies cluster expansion and streaming throughput.
"""

import sys
import os
import time
from pathlib import Path
from Bio import SeqIO
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

REPO_ROOT = Path("/Users/sergei/Projects/TOGA_MEME/recombination")
sys.path.insert(0, str(REPO_ROOT))

ALIGN_DIR = "/Users/sergei/Dropbox/Work/Collaborations/PittH5N1/2024-08-27/results/alignments"
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
LENS = [2277, 2271, 2148, 1701, 1494, 1407, 756, 690]

BENCH_DIR = REPO_ROOT / "benchmarks/flu_reassortment_benchmark"
RES_DIR = BENCH_DIR / "results"
RES_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 95)
print("  DYNAMIC MEDOID RECRUITMENT EVALUATION (7,054 COMPLETE 8-SEGMENT GENOMES)")
print("=" * 95)

# Step 1: Ingest all 8 segment alignments
t0 = time.perf_counter()
print("\n[*] Ingesting 8 segment alignments from PittH5N1 archive...")
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

# Build byte arrays per segment
print("\n[*] Converting alignments to contiguous byte arrays...")
seg_byte_arrays = {}
for s in SEGMENTS:
    seg_byte_arrays[s] = np.array([
        np.frombuffer(seg_indices[s][st].encode("ascii"), dtype=np.uint8)
        for st in common_strains
    ])

# Host categorization helper
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
np.random.seed(42)
host_buckets = {}
for st, h in strain_hosts.items():
    host_buckets.setdefault(h, []).append(st)

calib_strains = []
for h, strains in host_buckets.items():
    n_sample = max(2, int(400 * (len(strains) / N_total)))
    n_sample = min(n_sample, len(strains))
    calib_strains.extend(np.random.choice(strains, size=n_sample, replace=False).tolist())

remaining = [s for s in common_strains if s not in calib_strains]
if len(calib_strains) < 400:
    calib_strains.extend(np.random.choice(remaining, size=400 - len(calib_strains), replace=False).tolist())
calib_strains = sorted(calib_strains[:400])

strain_to_idx = {st: i for i, st in enumerate(common_strains)}
calib_indices = np.array([strain_to_idx[st] for st in calib_strains], dtype=int)

# Step 2: Discover Initial Reference Medoids on Diversity Core
print("\n[*] Discovering initial Regime A segment medoids via silhouette scoring...")
initial_medoids = {}
calibration_scores = {}

for s_name in SEGMENTS:
    all_bytes = seg_byte_arrays[s_name]
    calib_bytes = all_bytes[calib_indices]
    n_c = len(calib_indices)
    
    # Compute pairwise distance matrix on calibration core
    D_calib = np.zeros((n_c, n_c), dtype=float)
    for i in range(n_c):
        D_calib[i] = np.mean(calib_bytes != calib_bytes[i], axis=1)
        
    best_k = 2
    best_sil = -1.0
    best_labels = None
    for k in range(2, 6):
        c = AgglomerativeClustering(n_clusters=k, metric="precomputed", linkage="average")
        labs = c.fit_predict(D_calib)
        s = float(silhouette_score(D_calib, labs, metric="precomputed"))
        if s > best_sil:
            best_sil = s
            best_k = k
            best_labels = labs
            
    med_local_indices = []
    med_strains = []
    med_arrays = []
    for cl in range(best_k):
        members = np.where(best_labels == cl)[0]
        sub_D = D_calib[np.ix_(members, members)]
        best_local = members[int(np.argmin(sub_D.sum(axis=1)))]
        med_local_indices.append(best_local)
        med_strains.append(calib_strains[best_local])
        med_arrays.append(calib_bytes[best_local])
        
    # Non-conformity scores on calibration set
    med_mat = np.array(med_arrays) # shape (K, L_s)
    # dist from each calib sample to nearest medoid
    dists_to_meds = np.zeros((n_c, best_k), dtype=float)
    for c_id in range(best_k):
        dists_to_meds[:, c_id] = np.mean(calib_bytes != med_mat[c_id], axis=1)
    non_conf = np.sort(np.min(dists_to_meds, axis=1))
    
    initial_medoids[s_name] = {
        "strains": list(med_strains),
        "arrays": list(med_arrays), # list of (L_s,)
    }
    calibration_scores[s_name] = non_conf
    print(f"  - {s_name}: K = {best_k} medoids (Silhouette = {best_sil:.4f})")

# Step 3: Run Static Screening
print("\n[*] Evaluating Static Screening across 7,054 genomes...")
t_static = time.perf_counter()
static_p_vals = np.zeros((N_total, len(SEGMENTS)), dtype=float)
static_min_dists = np.zeros((N_total, len(SEGMENTS)), dtype=float)

for s_idx, s_name in enumerate(SEGMENTS):
    all_bytes = seg_byte_arrays[s_name]
    med_mat = np.array(initial_medoids[s_name]["arrays"])
    calib = calibration_scores[s_name]
    n_calib = len(calib)
    
    # Distance of all N_total to medoids
    d_to_meds = np.zeros((N_total, len(med_mat)), dtype=float)
    for m_id, m_arr in enumerate(med_mat):
        d_to_meds[:, m_id] = np.mean(all_bytes != m_arr, axis=1)
        
    min_d = np.min(d_to_meds, axis=1)
    # Conformal p-value
    counts_less = np.searchsorted(calib, min_d, side="left")
    p_vals = (1.0 + (n_calib - counts_less)) / (n_calib + 1.0)
    static_p_vals[:, s_idx] = p_vals
    static_min_dists[:, s_idx] = min_d

# Simes test
sorted_p_static = np.sort(static_p_vals, axis=1)
k_factors = 8.0 / np.arange(1, 9, dtype=float)
simes_static = np.min(sorted_p_static * k_factors[np.newaxis, :], axis=1)
static_novel_count = np.sum(simes_static < 0.05)
time_static = time.perf_counter() - t_static
print(f"[✓] Static Screen: {static_novel_count} novel genomes ({100.0 * static_novel_count / N_total:.2f}%) in {time_static:.4f} s")

# Step 4: Run Dynamic Medoid Recruitment
print("\n[*] Evaluating Dynamic Medoid Recruitment across 7,054 genomes...")
t_dyn = time.perf_counter()

# Deep copy of active medoids
dynamic_medoids = {
    s: {
        "strains": list(initial_medoids[s]["strains"]),
        "arrays": [arr.copy() for arr in initial_medoids[s]["arrays"]],
    }
    for s in SEGMENTS
}

alpha_recruit = 0.05 / 8.0  # Bonferroni segment threshold: 0.00625
min_recruit_dist = 0.005

recruited_events = []
dynamic_p_vals = np.zeros((N_total, len(SEGMENTS)), dtype=float)
dynamic_min_dists = np.zeros((N_total, len(SEGMENTS)), dtype=float)
dynamic_calls = []

# Streaming simulation over time / genome order
for i in range(N_total):
    strain_name = common_strains[i]
    q_p = np.zeros(len(SEGMENTS), dtype=float)
    q_d = np.zeros(len(SEGMENTS), dtype=float)
    genome_calls = []
    recruited_this = []
    
    for s_idx, s_name in enumerate(SEGMENTS):
        q_bytes = seg_byte_arrays[s_name][i]
        curr_meds = dynamic_medoids[s_name]["arrays"]
        calib = calibration_scores[s_name]
        n_calib = len(calib)
        
        # Distances to current active medoids
        dists = [float(np.mean(q_bytes != m)) for m in curr_meds]
        min_idx = int(np.argmin(dists))
        min_d = dists[min_idx]
        
        # Conformal p-value
        c_less = int(np.searchsorted(calib, min_d, side="left"))
        p_val = (1.0 + (n_calib - c_less)) / (n_calib + 1.0)
        
        q_p[s_idx] = p_val
        q_d[s_idx] = min_d
        
        # Dynamic recruitment check
        if p_val <= alpha_recruit and min_d >= min_recruit_dist:
            # Promote query to a new active medoid!
            dynamic_medoids[s_name]["arrays"].append(q_bytes)
            dynamic_medoids[s_name]["strains"].append(strain_name)
            new_cid = len(dynamic_medoids[s_name]["arrays"])
            
            recruited_events.append({
                "Step": i,
                "Taxon": strain_name,
                "Host": strain_hosts[strain_name],
                "Segment": s_name,
                "P_Value": p_val,
                "Min_Dist": min_d,
                "New_Cluster": f"C{new_cid}",
                "Founder_Medoid": strain_name
            })
            recruited_this.append(s_name)
            genome_calls.append(f"C{new_cid}*")
        else:
            genome_calls.append(f"C{min_idx + 1}")
            
    dynamic_p_vals[i] = q_p
    dynamic_min_dists[i] = q_d
    dynamic_calls.append("-".join(genome_calls))

# Simes test for dynamic run
sorted_p_dyn = np.sort(dynamic_p_vals, axis=1)
simes_dyn = np.min(sorted_p_dyn * k_factors[np.newaxis, :], axis=1)
dynamic_novel_count = np.sum(simes_dyn < 0.05)
time_dynamic = time.perf_counter() - t_dyn

print(f"[✓] Dynamic Screen: {dynamic_novel_count} novel genomes ({100.0 * dynamic_novel_count / N_total:.2f}%) in {time_dynamic:.4f} s")
print(f"    Recruited {len(recruited_events)} new exemplar medoids across all segments.")
print(f"    Alert reduction: {static_novel_count} -> {dynamic_novel_count} ({100.0 * (static_novel_count - dynamic_novel_count) / static_novel_count:.1f}% reduction in alarm fatigue)")

# Summarize recruited medoids
df_rec = pd.DataFrame(recruited_events)
df_rec.to_csv(RES_DIR / "dynamic_recruited_medoids.tsv", sep="\t", index=False)

print("\n[*] Recruited Medoids Breakdown by Segment:")
for s_name in SEGMENTS:
    sub = df_rec[df_rec["Segment"] == s_name] if len(df_rec) > 0 else []
    init_k = len(initial_medoids[s_name]["strains"])
    n_rec = len(sub)
    print(f"    - {s_name}: Initial {init_k} medoids -> Recruited {n_rec} -> Final {init_k + n_rec} medoids")

# Exemplar founder strains
if len(df_rec) > 0:
    print("\n[*] Top Recruited Founder Strains (Spanning Multiple Segments):")
    founder_counts = df_rec["Taxon"].value_counts().head(10)
    for st, count in founder_counts.items():
        h = strain_hosts[st]
        sub = df_rec[df_rec["Taxon"] == st]
        segs = ",".join(sub["Segment"].tolist())
        print(f"    - {st} [{h}]: Recruited for {count} segments ({segs})")

# Save dynamic typing table
df_dyn_out = pd.DataFrame({
    "Taxon": common_strains,
    "Host": [strain_hosts[st] for st in common_strains],
    "Constellation": dynamic_calls,
    "Omnibus_Simes_P": simes_dyn,
    "Is_Novel_Genome": simes_dyn < 0.05
})
for s_idx, s_name in enumerate(SEGMENTS):
    df_dyn_out[f"{s_name}_P"] = dynamic_p_vals[:, s_idx]
    df_dyn_out[f"{s_name}_MinDist"] = dynamic_min_dists[:, s_idx]

df_dyn_out.to_csv(RES_DIR / "dynamic_recruitment_audit.tsv", sep="\t", index=False)
print(f"\n[✓] Results successfully exported to {RES_DIR}")
