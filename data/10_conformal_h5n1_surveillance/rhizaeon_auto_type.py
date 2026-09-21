#!/usr/bin/env python3
"""
RhizAeon Autonomous Reassortment Typing Pipeline (rhizaeon_auto_type.py)

Implements the Two-Stage Active Surveillance Architecture:
- Stage 1 (Regime A): Diversity Core Sampling and Unsupervised Segment Medoid Discovery.
- Stage 2 (Calibration): Conformal Prediction Calibration for Non-Arbitrary Novelty Detection.
- Stage 3 (Regime B): Vectorized Linear-Time Whole-Genome Constellation Typing.
- Stage 4 (Outlier Audit): Conformal p-value assignment and Novelty Buffer reporting.
"""

import sys
import os
import time
import argparse
from pathlib import Path
from Bio import SeqIO
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

REPO_ROOT = Path("/Users/sergei/Projects/TOGA_MEME/recombination")
sys.path.insert(0, str(REPO_ROOT))

from rhizaeon.tensor import encode_alignment_matrix, PrefixDistanceEngine

DEFAULT_SEGMENTS = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]
DEFAULT_LENS = [2277, 2271, 2148, 1701, 1494, 1407, 756, 690]

def select_diversity_core(engine, taxa, n_core=120):
    """
    Minimax Furthest-Point Sampling (FPS) in Whole-Genome Manifold Space.
    Guarantees maximal geometric dispersion across sequence space.
    """
    N = len(taxa)
    if N <= n_core:
        return list(range(N))
        
    D_full = engine.query_distance_matrix(0, engine.L)
    
    # 1. Start with the taxon with maximum variance in pairwise distance
    first_idx = int(np.argmax(np.var(D_full, axis=1)))
    selected = [first_idx]
    
    # 2. Track minimum distance from each candidate to the selected set
    min_dists = D_full[first_idx].copy()
    
    # 3. Iteratively pick the point that maximizes the minimum distance
    for _ in range(1, n_core):
        next_idx = int(np.argmax(min_dists))
        selected.append(next_idx)
        min_dists = np.minimum(min_dists, D_full[next_idx])
        
    return sorted(selected)

def discover_segment_medoids(engine, core_indices, taxa, segment_bounds, k_range=(2, 5)):
    """
    Regime A: Discovers optimal cluster count K_k and exemplar sequence medoids
    independently for each segment via silhouette optimization.
    """
    n_core = len(core_indices)
    discovered_centroids = {}
    
    for s_idx, (start, end) in enumerate(segment_bounds):
        # Extract segment distance matrix for the core in O(1)
        D_seg_full = engine.query_distance_matrix(start, end)
        D_core = D_seg_full[np.ix_(core_indices, core_indices)]
        
        best_k = 2
        best_sil = -1.0
        best_labels = np.zeros(n_core, dtype=int)
        
        for k in range(k_range[0], min(k_range[1] + 1, n_core)):
            clusterer = AgglomerativeClustering(
                n_clusters=k, metric="precomputed", linkage="average"
            )
            labels = clusterer.fit_predict(D_core)
            sil = float(silhouette_score(D_core, labels, metric="precomputed"))
            if sil > best_sil:
                best_sil = sil
                best_k = k
                best_labels = labels
                
        # Discover exemplar medoids and non-conformity calibration scores
        medoid_taxa = []
        medoid_core_indices = []
        cluster_radii = []
        non_conformity_scores = []
        
        for c in range(best_k):
            members = np.where(best_labels == c)[0]
            sub_D = D_core[np.ix_(members, members)]
            # Medoid minimizes sum of distances to members
            best_local_idx = int(np.argmin(sub_D.sum(axis=1)))
            medoid_core_idx = members[best_local_idx]
            global_medoid_idx = core_indices[medoid_core_idx]
            
            medoid_taxa.append(taxa[global_medoid_idx])
            medoid_core_indices.append(medoid_core_idx)
            
            # Radii and member distances
            dists_to_medoid = sub_D[best_local_idx]
            cluster_radii.append(float(np.max(dists_to_medoid)))
            non_conformity_scores.extend(dists_to_medoid.tolist())
            
        # Conformal calibration: 99th percentile of intra-cluster distance
        conformal_threshold = float(np.percentile(non_conformity_scores, 99.0))
        
        discovered_centroids[s_idx] = {
            "optimal_k": best_k,
            "silhouette": best_sil,
            "medoids": medoid_taxa,
            "medoid_core_indices": medoid_core_indices,
            "global_medoid_indices": [core_indices[i] for i in medoid_core_indices],
            "radii": cluster_radii,
            "conformal_threshold": conformal_threshold,
            "calibration_scores": np.array(non_conformity_scores)
        }
        
    return discovered_centroids

def stream_regime_b_typing(engine, taxa, segment_bounds, discovered_centroids):
    """
    Regime B: Vectorized Linear-Time Constellation Calling and Conformal Outlier Detection.
    Runs in O(N * K_max * L) memory and time.
    """
    N = len(taxa)
    records = []
    
    for s_idx, (start, end) in enumerate(segment_bounds):
        info = discovered_centroids[s_idx]
        ref_indices = info["global_medoid_indices"]
        K = len(ref_indices)
        
        # Query distance from all N taxa to the K medoids in O(1) prefix operations
        D_seg_all = engine.query_distance_matrix(start, end)
        D_to_refs = D_seg_all[:, ref_indices] # Shape (N, K)
        
        # Closest centroid
        best_cluster = np.argmin(D_to_refs, axis=1)
        min_dists = np.min(D_to_refs, axis=1)
        
        # Conformal empirical p-value: fraction of calibration distances >= observed
        calib = info["calibration_scores"]
        # Vectorized p-value estimation: P(D_calib >= d_query)
        # Using searchsorted on sorted calibration scores
        calib_sorted = np.sort(calib)
        p_values = 1.0 - (np.searchsorted(calib_sorted, min_dists, side="left") / len(calib_sorted))
        
        # Store segment calls
        info["calls"] = best_cluster
        info["distances"] = min_dists
        info["p_values"] = p_values
        
    # Compile whole-genome constellation records
    for i in range(N):
        calls = [f"C{discovered_centroids[s]['calls'][i] + 1}" for s in range(len(segment_bounds))]
        min_p = min(discovered_centroids[s]["p_values"][i] for s in range(len(segment_bounds)))
        outlier_segs = [
            s for s in range(len(segment_bounds))
            if discovered_centroids[s]["p_values"][i] < 0.01
        ]
        
        records.append({
            "Taxon": taxa[i],
            "Constellation": "-".join(calls),
            "Min_Conformal_P": float(min_p),
            "Is_Novel_Outlier": len(outlier_segs) > 0,
            "Outlier_Segments": ",".join(str(s+1) for s in outlier_segs) if outlier_segs else "None"
        })
        
    return pd.DataFrame(records)

def stream_dynamic_regime_b_typing(engine, taxa, segment_bounds, discovered_centroids, alpha_recruit=0.01, min_recruit_dist=0.005):
    """
    Regime B with Dynamic Medoid Recruitment (Online Adaptive Regime A).
    Promotes statistically verified novel segment alleles into the active
    exemplar reference medoid library during streaming inference.
    """
    N = len(taxa)
    n_segments = len(segment_bounds)
    
    # Pre-extract segment distance matrices in O(1)
    seg_dists = [
        engine.query_distance_matrix(start, end)
        for (start, end) in segment_bounds
    ]
    
    # Active reference medoid indices per segment
    active_indices = [
        list(discovered_centroids[s]["global_medoid_indices"])
        for s in range(n_segments)
    ]
    calibrations = [
        np.sort(discovered_centroids[s]["calibration_scores"])
        for s in range(n_segments)
    ]
    
    recruited_events = []
    records = []
    
    for i in range(N):
        calls = []
        p_values = []
        outlier_segs = []
        recruited_segs = []
        
        for s in range(n_segments):
            D_s = seg_dists[s]
            curr_meds = active_indices[s]
            dists = D_s[i, curr_meds]
            min_idx = int(np.argmin(dists))
            min_d = float(dists[min_idx])
            
            calib = calibrations[s]
            n_calib = len(calib)
            c_less = int(np.searchsorted(calib, min_d, side="left"))
            p_val = (1.0 + (n_calib - c_less)) / (n_calib + 1.0)
            p_values.append(p_val)
            
            if p_val <= alpha_recruit and min_d >= min_recruit_dist:
                # Promote to new active medoid
                active_indices[s].append(i)
                new_cid = len(active_indices[s])
                recruited_events.append({
                    "Step": i,
                    "Taxon": taxa[i],
                    "Segment_Idx": s + 1,
                    "Segment_Name": DEFAULT_SEGMENTS[s],
                    "P_Value": p_val,
                    "Min_Dist": min_d,
                    "Cluster_ID": f"C{new_cid}"
                })
                recruited_segs.append(str(s + 1))
                calls.append(f"C{new_cid}*")
                outlier_segs.append(s)
            else:
                calls.append(f"C{min_idx + 1}")
                if p_val <= alpha_recruit:
                    outlier_segs.append(s)
                    
        records.append({
            "Taxon": taxa[i],
            "Constellation": "-".join(calls),
            "Min_Conformal_P": float(min(p_values)),
            "Is_Novel_Outlier": len(outlier_segs) > 0,
            "Outlier_Segments": ",".join(str(s+1) for s in outlier_segs) if outlier_segs else "None",
            "Recruited_Segments": ",".join(recruited_segs) if recruited_segs else "None"
        })
        
    return pd.DataFrame(records), pd.DataFrame(recruited_events)

def main():
    parser = argparse.ArgumentParser(description="RhizAeon Autonomous Two-Stage Influenza A Typing Pipeline")
    parser.add_argument("--fasta", type=str, required=True, help="Path to concatenated pseudomolecule FASTA")
    parser.add_argument("--core_size", type=int, default=120, help="Size of Regime A diversity core")
    parser.add_argument("--dynamic_recruitment", action="store_true", help="Enable online adaptive dynamic medoid recruitment")
    parser.add_argument("--out_catalog", type=str, default="autonomous_reassortment_catalog.tsv", help="Output catalog TSV")
    parser.add_argument("--out_summary", type=str, default="autonomous_reference_manifest.csv", help="Discovered references CSV")
    parser.add_argument("--out_recruited", type=str, default="autonomous_recruited_medoids.tsv", help="Recruited medoids TSV")
    args = parser.parse_args()
    
    print("=" * 95)
    print("  RHIZAEON AUTONOMOUS TWO-STAGE REASSORTMENT TYPING ENGINE")
    print("=" * 95)
    
    t0 = time.perf_counter()
    print(f"\n[*] Encoding alignment: {args.fasta}...")
    seq_mat, taxa, L = encode_alignment_matrix(args.fasta)
    engine = PrefixDistanceEngine(seq_mat, codon_aligned=False)
    print(f"[✓] Alignment Encoded: {len(taxa):,} taxa, {L:,} nucleotides in {time.perf_counter()-t0:.2f} s")
    
    # Segment boundaries
    cum = 0
    bounds = []
    for l in DEFAULT_LENS:
        bounds.append((cum, cum + l))
        cum += l
        
    # Stage 1: Diversity core selection
    print(f"\n[*] Stage 1: Selecting diversity core (N_0 = {args.core_size}) via Furthest-Point Sampling...")
    t0 = time.perf_counter()
    core_idx = select_diversity_core(engine, taxa, n_core=args.core_size)
    print(f"[✓] Diversity core selected in {time.perf_counter()-t0:.2f} s")
    
    # Stage 2: Regime A Medoid Discovery
    print("\n[*] Stage 2: Discovering segment-specific reference medoids and calibrating conformal scores...")
    t0 = time.perf_counter()
    centroids = discover_segment_medoids(engine, core_idx, taxa, bounds, k_range=(2, 6))
    print(f"[✓] Reference centroids discovered in {time.perf_counter()-t0:.2f} s")
    
    manifest_rows = []
    for s_idx, name in enumerate(DEFAULT_SEGMENTS):
        info = centroids[s_idx]
        print(f"    • {name} (span {bounds[s_idx][0]}-{bounds[s_idx][1]} nt): Optimal K = {info['optimal_k']} (Silhouette: {info['silhouette']:.3f}, Conformal Threshold: {info['conformal_threshold']:.4f})")
        for m_id, m_taxon in enumerate(info["medoids"]):
            clean_name = m_taxon.split("|")[2] if "|" in m_taxon else m_taxon
            print(f"        - Centroid {m_id+1}: {clean_name} (Radius: {info['radii'][m_id]:.4f})")
            manifest_rows.append({
                "Segment": name,
                "Cluster_ID": m_id + 1,
                "Medoid_Taxon": m_taxon,
                "Cluster_Radius": info["radii"][m_id],
                "Conformal_Threshold": info["conformal_threshold"],
                "Silhouette": info["silhouette"]
            })
            
    pd.DataFrame(manifest_rows).to_csv(args.out_summary, index=False)
    print(f"[✓] Discovered reference manifest saved to {args.out_summary}")
    
    # Stage 3: Regime B Linear Typing
    if args.dynamic_recruitment:
        print(f"\n[*] Stage 3: Streaming Regime B with DYNAMIC MEDOID RECRUITMENT across {len(taxa):,} taxa...")
        t0 = time.perf_counter()
        df_results, df_recruited = stream_dynamic_regime_b_typing(engine, taxa, bounds, centroids)
        t_typing = time.perf_counter() - t0
        print(f"[✓] Completed whole-genome dynamic typing in {t_typing:.3f} s ({t_typing/len(taxa)*1000:.3f} ms/genome)")
        print(f"[✓] Dynamically recruited {len(df_recruited)} new exemplar medoids during streaming.")
        df_recruited.to_csv(args.out_recruited, sep="\t", index=False)
        print(f"[✓] Recruited medoid audit manifest saved to {args.out_recruited}")
    else:
        print(f"\n[*] Stage 3: Streaming Regime B Static Whole-Genome Typing across {len(taxa):,} taxa...")
        t0 = time.perf_counter()
        df_results = stream_regime_b_typing(engine, taxa, bounds, centroids)
        t_typing = time.perf_counter() - t0
        print(f"[✓] Completed whole-genome static typing in {t_typing:.3f} s ({t_typing/len(taxa)*1000:.3f} ms/genome)")
        
    df_results.to_csv(args.out_catalog, sep="\t", index=False)
    print(f"[✓] Full classification catalog saved to {args.out_catalog}")
    
    # Summary of findings
    n_distinct = df_results["Constellation"].nunique()
    n_outliers = int(df_results["Is_Novel_Outlier"].sum())
    print(f"\nSummary of Epizootic Classification:")
    print(f"  • Total Genomes Classified: {len(taxa):,}")
    print(f"  • Distinct Reassortant Constellations: {n_distinct}")
    print(f"  • Statistically Verified Novel Outliers (Conformal p < 0.01): {n_outliers} ({n_outliers/len(taxa)*100:.2f}%)")
    print(f"\nTop 5 Discovered Constellations:")
    print(df_results["Constellation"].value_counts().head(5).to_string())

if __name__ == "__main__":
    main()

