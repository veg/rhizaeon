#!/usr/bin/env python3
"""
Direct De Novo Reassortment Discovery Screen across 120 Complete Influenza A Genomes using Modern RhizAeon.
Analyzes 12,744 nt concatenated pseudomolecules across:
- Human zoonotic cases (Texas/37/2024, Michigan/90/2024)
- Bovine (dairy cattle B3.13 across TX, KS, ID, MI, SD, NM)
- Feline (domestic cats)
- Terrestrial carnivores (red foxes, striped skunks, black bears)
- Marine mammals (harbor seals, South American sea lions)
- Avian raptors (bald eagles, peregrine falcons, hawks)
- Seabirds & shorebirds (terns, gulls)
- Wild waterfowl (teals, wigeons, mallards, brants, geese)
- Domestic poultry (chickens, turkeys)
- Eurasian reference lineages
"""

import sys
import os
import time
from pathlib import Path
from Bio import SeqIO
import pandas as pd
import numpy as np

REPO_ROOT = Path("/Users/sergei/Projects/TOGA_MEME/recombination")
sys.path.insert(0, str(REPO_ROOT))

from rhizaeon.tensor import encode_alignment_matrix, PrefixDistanceEngine
from rhizaeon.fda import run_fda_recombination_screen, run_recursive_partition_fda_screen

ALIGN_DIR = "/Users/sergei/Dropbox/Work/Collaborations/PittH5N1/2024-08-27/results/alignments"
SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
LENS = [2277, 2271, 2148, 1701, 1494, 1407, 756, 690]
TOTAL_L = sum(LENS) # 12744 nt

BENCH_DIR = REPO_ROOT / "benchmarks/flu_reassortment_benchmark"
DATA_DIR = BENCH_DIR / "data"
RES_DIR = BENCH_DIR / "results"
DATA_DIR.mkdir(parents=True, exist_ok=True)
RES_DIR.mkdir(parents=True, exist_ok=True)

print("=" * 95)
print("  DIRECT DE NOVO REASSORTMENT DISCOVERY SCREEN WITH RHIZAEON (12,744 nt PSEUDOMOLECULES)")
print("=" * 95)

# Step 1: Index all 8 segments
print("\n[*] Indexing 8 segment alignments from PittH5N1 archive...")
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
    print(f"    • Indexed {len(d)} records for {s} (length: {LENS[len(seg_indices)-1]} nt)")

common_strains = sorted(list(set.intersection(*[set(d.keys()) for d in seg_indices.values()])))
print(f"[✓] Found {len(common_strains):,} strains with complete, uninterrupted 8-segment sequences in {time.perf_counter()-t0:.2f} s")

# Step 2: Stratified selection of 120 representative genomes
def categorize_strain(st):
    s = st.lower()
    if "texas/37" in s or "michigan/90" in s:
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
    if "dairy_cow" in s or "cattle" in s:
        return "Bovine"
    if "eagle" in s or "falcon" in s or "hawk" in s:
        return "Raptor"
    if "gull" in s or "tern" in s or "sandpiper" in s:
        return "Seabird_Shorebird"
    if "teal" in s or "wigeon" in s or "mallard" in s or "brant" in s or "duck" in s or "goose" in s:
        return "Waterfowl"
    if "chicken" in s or "turkey" in s:
        return "Poultry"
    return "Other"

# Target counts per category
target_quotas = {
    "Human": 5,
    "Bear": 5,
    "Marine_Mammal": 12,
    "Cat": 10,
    "Skunk": 10,
    "Fox": 10,
    "Bovine": 15,
    "Raptor": 12,
    "Seabird_Shorebird": 12,
    "Waterfowl": 22,
    "Poultry": 12
}

selected_strains = []
# Ensure key anchors are explicitly included:
anchors = [
    "A/goose/Virginia/W22-773/2022",
    "A/blue-winged_teal/Texas/UGAI22-2961/2022",
    "A/American_Wigeon/North_Carolina/USDA-001330-028/2022",
    "A/dairy_cow/USA/24_019691-004/2024",
    "A/Texas/37/2024",
    "A/Michigan/90/2024"
]
for a in anchors:
    if a in common_strains and a not in selected_strains:
        selected_strains.append(a)

st_by_cat = {c: [] for c in target_quotas}
for st in common_strains:
    c = categorize_strain(st)
    if c in st_by_cat:
        st_by_cat[c].append(st)

np.random.seed(42)
for c, q in target_quotas.items():
    cand = [s for s in st_by_cat[c] if s not in selected_strains]
    n_pick = min(len(cand), q - sum(1 for s in selected_strains if categorize_strain(s) == c))
    if n_pick > 0:
        picked = np.random.choice(cand, size=n_pick, replace=False).tolist()
        selected_strains.extend(picked)

print(f"\n[✓] Assembled discovery cohort: {len(selected_strains)} complete 8-segment genomes across 11 host clades.")

# Step 3: Build concatenated alignment (12,744 nt)
concat_fa = DATA_DIR / "concatenated_120_genomes.fasta"
print(f"[*] Constructing 12,744 nt pseudomolecules and writing to {concat_fa}...")
with open(concat_fa, "w") as fp:
    for st in selected_strains:
        full_seq = "".join(seg_indices[s][st] for s in SEGMENTS)
        clean_name = st.replace("/", "_").replace("-", "_").replace(" ", "_")
        fp.write(f">{clean_name}\n{full_seq}\n")

# Step 4: Run RhizAeon Prefix Distance Engine
t0 = time.perf_counter()
seq_mat, taxa, L = encode_alignment_matrix(str(concat_fa))
t_enc = (time.perf_counter() - t0) * 1000
N = len(taxa)

t0 = time.perf_counter()
engine = PrefixDistanceEngine(seq_mat, codon_aligned=False)
t_pref = (time.perf_counter() - t0) * 1000
print(f"[✓] Alignment Encoded (N={N}, L={L:,} nt) in {t_enc:.1f} ms | Prefix Distance Tensor built in {t_pref:.1f} ms")

# Step 5: Segment boundaries
cum = 0
bounds = [0]
for l in LENS:
    cum += l
    bounds.append(cum)

# Find reference indices
ref_eur = next(i for i, t in enumerate(selected_strains) if "W22-773" in t)
ref_am1 = next(i for i, t in enumerate(selected_strains) if "UGAI22-2961" in t)
ref_am2 = next(i for i, t in enumerate(selected_strains) if "USDA-001330" in t)
ref_cow = next(i for i, t in enumerate(selected_strains) if "24_019691-004" in t)

print(f"\n[*] Reference Anchors in Manifold:")
print(f"    • Ancestral Eurasian 2.3.4.4b Anchor: {selected_strains[ref_eur]}")
print(f"    • North American Wild Bird LPAI-1:   {selected_strains[ref_am1]}")
print(f"    • North American Wild Bird LPAI-2:   {selected_strains[ref_am2]}")
print(f"    • Dairy Cattle B3.13 Anchor:         {selected_strains[ref_cow]}")

# Step 6: O(1) Segment Distance Queries and Constellation Calling
print("\n[*] Calling 8-segment constellation vectors via O(1) Prefix Distance Tensors...")
t0 = time.perf_counter()

constellation_records = []
for i, st in enumerate(selected_strains):
    host = categorize_strain(st)
    seg_types = []
    seg_dists_eur = []
    seg_dists_am1 = []
    seg_dists_am2 = []
    
    for s_idx, s in enumerate(SEGMENTS):
        b_start = bounds[s_idx]
        b_end = bounds[s_idx + 1]
        D = engine.query_distance_matrix(b_start, b_end)
        
        d_e = D[i, ref_eur]
        d_a1 = D[i, ref_am1]
        d_a2 = D[i, ref_am2]
        
        seg_dists_eur.append(d_e)
        seg_dists_am1.append(d_a1)
        seg_dists_am2.append(d_a2)
        
        # Classification rule:
        # If distance to Eurasian reference is very low (<0.025), classified as Eurasian (EA)
        # Else if closer to North American LPAI 1 than Eurasian, classified as AM1
        # Else if closer to North American LPAI 2, classified as AM2
        min_d = min(d_e, d_a1, d_a2)
        if min_d == d_e and d_e < 0.030:
            assigned = "EA"
        elif min_d == d_a1:
            assigned = "AM1"
        elif min_d == d_a2:
            assigned = "AM2"
        else:
            assigned = "EA"
        seg_types.append(assigned)
        
    constellation_code = "-".join(seg_types)
    constellation_records.append({
        "Strain": st,
        "Host_Group": host,
        "Constellation": constellation_code,
        "PB2": seg_types[0],
        "PB1": seg_types[1],
        "PA": seg_types[2],
        "HA": seg_types[3],
        "NP": seg_types[4],
        "NA": seg_types[5],
        "MP": seg_types[6],
        "NS": seg_types[7]
    })

t_typing = (time.perf_counter() - t0) * 1000
df_const = pd.DataFrame(constellation_records)
df_const.to_csv(RES_DIR / "de_novo_reassortment_constellations.csv", index=False)
print(f"[✓] Constellation typing finished for all {N} genomes in {t_typing:.1f} ms ({t_typing/N:.2f} ms/genome)!")

# Step 7: Unsupervised Continuous Trajectory Screen across concatenated genomes
print("\n[*] Executing Unsupervised Continuous Trajectory Screen across 12,744 nt...")
t0 = time.perf_counter()
_, bps_traj = run_fda_recombination_screen(
    engine, taxa,
    k=8,
    bin_size=400,
    step=40,
    min_kinetic_z=2.2,
    min_pir=0.04,
    polish_ml=True
)
t_screen = (time.perf_counter() - t0) * 1000
print(f"[✓] Continuous Sequence Manifold Screen complete in {t_screen:.1f} ms")
print(f"    Detected {len(bps_traj)} unsupervised changepoints across chromosome.")

bp_records = []
known_junctions = [2277, 4548, 6696, 8397, 9891, 11298, 12054]
junction_names = ["PB2/PB1", "PB1/PA", "PA/HA", "HA/NP", "NP/NA", "NA/MP", "MP/NS"]

for b in bps_traj:
    # Match closest segment junction
    best_j_idx = int(np.argmin([abs(b.breakpoint_nt - j) for j in known_junctions]))
    j_nt = known_junctions[best_j_idx]
    j_name = junction_names[best_j_idx]
    delta = abs(b.breakpoint_nt - j_nt)
    is_boundary = delta <= 150
    bp_records.append({
        "Inferred_BP": b.breakpoint_nt,
        "CI_Left": b.ci_left,
        "CI_Right": b.ci_right,
        "Plateau_Width": b.plateau_width,
        "Closest_Junction": j_name,
        "Junction_Coordinate": j_nt,
        "Delta_nt": delta,
        "Is_Reassortment_Boundary": is_boundary,
        "Recombinant_Taxon": b.recombinant_taxon,
        "Parent_1": b.parent_1,
        "Parent_2": b.parent_2,
        "Kinetic_Z": round(b.kinetic_z, 2),
        "dLL": round(b.log_likelihood_gain, 1) if b.log_likelihood_gain else 0.0
    })

df_bps = pd.DataFrame(bp_records)
df_bps.to_csv(RES_DIR / "unsupervised_breakpoints.csv", index=False)

# Step 8: Constellation Diversity Analysis and Novel Animal Reassortants
print("\n" + "=" * 95)
print("  REASSORTMENT CONSTELLATION DIVERSITY ACROSS HOST CLADES")
print("=" * 95)

const_counts = df_const["Constellation"].value_counts()
print(f"Discovered {len(const_counts)} distinct reassortant genotypes across {N} genomes:")
for c_code, count in const_counts.head(10).items():
    hosts = df_const[df_const["Constellation"] == c_code]["Host_Group"].unique().tolist()
    print(f"  • Genotype [{c_code}]: {count:2d} isolates | Hosts: {', '.join(hosts)}")

# Host breakdown
print("\nHost-Stratified Genotype Summary:")
host_summary = []
for h, grp in df_const.groupby("Host_Group"):
    top_c = grp["Constellation"].mode()[0]
    n_top = sum(grp["Constellation"] == top_c)
    n_distinct = grp["Constellation"].nunique()
    host_summary.append({
        "Host_Group": h,
        "Total_Genomes": len(grp),
        "Distinct_Genotypes": n_distinct,
        "Dominant_Genotype": top_c,
        "Dominant_Prevalence": f"{n_top}/{len(grp)} ({n_top/len(grp)*100:.1f}%)"
    })
    print(f"  {h:<18s}: {len(grp):2d} genomes | {n_distinct:2d} distinct genotypes | Dominant: [{top_c}] ({n_top}/{len(grp)})")

df_host = pd.DataFrame(host_summary)
df_host.to_csv(RES_DIR / "host_genotype_summary.csv", index=False)

# Step 9: Save Novel Reassortant Highlights
# Find non-dominant mammalian reassortants
bovine_dominant = df_const[df_const["Host_Group"] == "Bovine"]["Constellation"].mode()[0]
print(f"\n[*] Dominant Bovine B3.13 Genotype Constellation: [{bovine_dominant}]")

divergent_mammals = df_const[(df_const["Host_Group"].isin(["Bear", "Marine_Mammal", "Cat", "Skunk", "Fox"])) & (df_const["Constellation"] != bovine_dominant)]
print(f"[!] Discovered {len(divergent_mammals)} mammalian isolates carrying NON-B3.13 reassortant genotypes:")
for _, r in divergent_mammals.head(15).iterrows():
    print(f"    • {r['Strain']:<55s} | Host: {r['Host_Group']:<14s} | Constellation: [{r['Constellation']}]")

divergent_mammals.to_csv(RES_DIR / "novel_mammalian_reassortants.csv", index=False)

print("\n" + "=" * 95)
print(f"Discovery screen complete. All datasets, prefix distance metrics, and CSV logs archived in:")
print(f"  • {RES_DIR}")
print("=" * 95)
