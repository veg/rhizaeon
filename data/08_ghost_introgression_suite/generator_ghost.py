"""
benchmarks/multi_taxon_simulations/generator_multi.py
=====================================================
Multi-Taxon Synthetic Alignment Generators for Reticulate Evolution Benchmarks.

Generates parametric alignments across:
  Option A: Taxonomic Scaling (N in {5, 10, 20, 50, 100})
  Option B: Ghost Donors / Incomplete Sampling (delta in {0%, 1%, 2.5%, 5%, 10%, 20%})
  Option C: Circulating Recombinant Forms / Clade Expansion (M in {1, 2, 4, 8, 12})
"""

import os
import sys
import numpy as np
import pyvolve
from typing import Dict, List, Tuple, Optional

NT_CHARS = ["A", "C", "G", "T", "-"]
CHAR_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "-": 4, "N": 4}

def seqs_to_matrix(taxa: List[str], seqs: Dict[str, str], L: int) -> np.ndarray:
    mat = np.zeros((len(taxa), L), dtype=np.int8)
    for i, t in enumerate(taxa):
        s = seqs[t].upper()
        for j in range(min(L, len(s))):
            mat[i, j] = CHAR_MAP.get(s[j], 0)
    return mat

# -----------------------------------------------------------------------------
# Option A: Taxonomic Scaling (N in {5, 10, 20, 50, 100})
# -----------------------------------------------------------------------------

def generate_scaling_alignment(
    N: int = 20,
    divergence: float = 0.05,
    tract_len: int = 500,
    tract_start: int = 1000,
    L: int = 3000,
    seed: int = 42,
    is_null: bool = False
) -> Tuple[np.ndarray, List[str], int, Dict]:
    """
    Generates an N-taxon alignment structured into two divergent clades (Clade A and Clade B).
    If not is_null, taxon 'Rec' is a recombinant:
      Flanks [0, tract_start) and [tract_start + tract_len, L) belong to Clade A.
      Cassette [tract_start, tract_start + tract_len) belongs to Clade B.
    If is_null, all N taxa evolve strictly clonally under the same tree with zero recombination.
    """
    np.random.seed(seed)
    n_a = N // 2
    n_b = N - n_a

    # Sub-clade intra-divergence: ~0.2 * divergence
    d_within = max(0.005, divergence * 0.20)
    d_between = divergence

    a_taxa = [f"A{i}" for i in range(n_a)]
    b_taxa = [f"B{i}" for i in range(n_b)]

    if is_null:
        # Strictly clonal tree
        c1_str = ",".join(f"{t}:{d_within:.5f}" for t in a_taxa)
        c2_str = ",".join(f"{t}:{d_within:.5f}" for t in b_taxa)
        tree_str = f"(({c1_str}):{d_between/2:.5f}, ({c2_str}):{d_between/2:.5f});"
        tree = pyvolve.read_tree(tree=tree_str)

        model = pyvolve.Model("nucleotide", {"parameters": {"kappa": 2.5}})
        p = pyvolve.Partition(models=model, size=L)
        ev = pyvolve.Evolver(partitions=p, tree=tree)
        ev(ratefile=None, infofile=None, seqfile=None, seed=seed)
        seqs = ev.get_sequences()
        taxa = sorted(list(seqs.keys()))
        mat = seqs_to_matrix(taxa, seqs, L)
        meta = {"is_null": True, "N": N, "divergence": divergence, "true_bps": []}
        return mat, taxa, L, meta

    else:
        # Recombinant taxon 'Rec' switches parentage
        # Left/Right Flanks: Rec is in Clade A (sister to A0)
        c1_left = f"(Rec:{d_within:.5f}, A0:{d_within:.5f}):{d_within:.5f}," + ",".join(f"{t}:{d_within:.5f}" for t in a_taxa[1:])
        c2_left = ",".join(f"{t}:{d_within:.5f}" for t in b_taxa)
        tree_flank_str = f"(({c1_left}):{d_between/2:.5f}, ({c2_left}):{d_between/2:.5f});"
        tree_flank = pyvolve.read_tree(tree=tree_flank_str)

        # Cassette: Rec is in Clade B (sister to B0)
        c1_cass = ",".join(f"{t}:{d_within:.5f}" for t in a_taxa)
        c2_cass = f"(Rec:{d_within:.5f}, B0:{d_within:.5f}):{d_within:.5f}," + ",".join(f"{t}:{d_within:.5f}" for t in b_taxa[1:])
        tree_cass_str = f"(({c1_cass}):{d_between/2:.5f}, ({c2_cass}):{d_between/2:.5f});"
        tree_cass = pyvolve.read_tree(tree=tree_cass_str)

        model = pyvolve.Model("nucleotide", {"parameters": {"kappa": 2.5}})
        p1 = pyvolve.Partition(models=model, size=tract_start)
        p2 = pyvolve.Partition(models=model, size=tract_len)
        p3 = pyvolve.Partition(models=model, size=L - tract_start - tract_len)

        e1 = pyvolve.Evolver(partitions=p1, tree=tree_flank)
        e1(ratefile=None, infofile=None, seqfile=None, seed=seed)
        e2 = pyvolve.Evolver(partitions=p2, tree=tree_cass)
        e2(ratefile=None, infofile=None, seqfile=None, seed=seed + 1)
        e3 = pyvolve.Evolver(partitions=p3, tree=tree_flank)
        e3(ratefile=None, infofile=None, seqfile=None, seed=seed + 2)

        s1, s2, s3 = e1.get_sequences(), e2.get_sequences(), e3.get_sequences()
        taxa = sorted(list(s1.keys()))
        seqs = {t: s1[t] + s2[t] + s3[t] for t in taxa}
        mat = seqs_to_matrix(taxa, seqs, L)
        meta = {
            "is_null": False,
            "N": N,
            "divergence": divergence,
            "recombinant": "Rec",
            "parent_left": "A0",
            "parent_right": "B0",
            "true_bps": [tract_start, tract_start + tract_len],
            "tract_len": tract_len
        }
        return mat, taxa, L, meta


# -----------------------------------------------------------------------------
# Option B: Ghost Donors / Incomplete Lineage Sampling
# -----------------------------------------------------------------------------

def generate_ghost_alignment(
    divergence: float = 0.05,
    donor_divergence: float = 0.05,
    is_complete_ghost: bool = False,
    tract_len: int = 600,
    tract_start: int = 1000,
    L: int = 3000,
    seed: int = 42
) -> Tuple[np.ndarray, List[str], int, Dict]:
    """
    Evaluates recombination where the donor lineage (Parent 2) is either:
      - Present but diverged by donor_divergence from the true cassette donor.
      - A complete 'ghost' lineage (is_complete_ghost=True), where Clade B is entirely unrepresented.
    Total taxa: 12 taxa (or 8 if complete ghost).
    """
    np.random.seed(seed)
    d = divergence / 2.0
    d_donor = donor_divergence

    # Clade A: A0, A1, A2, A3
    # Clade B (Donor Clade): B0 (true donor), B1, B2, B3
    # Clade C (Outgroup): C0, C1, C2, C3
    taxa_a = ["A0", "A1", "A2", "A3"]
    taxa_b = ["B0", "B1", "B2", "B3"]
    taxa_c = ["C0", "C1", "C2", "C3"]

    c_a_str = ",".join(f"{t}:0.0100" for t in taxa_a)
    c_b_str = ",".join(f"{t}:0.0100" for t in taxa_b)
    c_c_str = ",".join(f"{t}:0.0100" for t in taxa_c)

    # Flanks: Rec is sister to A0
    tree_flank_str = f"(((Rec:0.0100, {c_a_str}):{d:.5f}, ({c_b_str}):{d:.5f}):{d:.5f}, ({c_c_str}):{2*d:.5f});"
    # Cassette: Rec is derived from B0, but B0 may have drifted by d_donor
    tree_cass_str = f"((({c_a_str}):{d:.5f}, (Rec:{d_donor:.5f}, {c_b_str}):{d:.5f}):{d:.5f}, ({c_c_str}):{2*d:.5f});"

    model = pyvolve.Model("nucleotide", {"parameters": {"kappa": 2.5}})
    p1 = pyvolve.Partition(models=model, size=tract_start)
    p2 = pyvolve.Partition(models=model, size=tract_len)
    p3 = pyvolve.Partition(models=model, size=L - tract_start - tract_len)

    e1 = pyvolve.Evolver(partitions=p1, tree=pyvolve.read_tree(tree=tree_flank_str))
    e1(ratefile=None, infofile=None, seqfile=None, seed=seed)
    e2 = pyvolve.Evolver(partitions=p2, tree=pyvolve.read_tree(tree=tree_cass_str))
    e2(ratefile=None, infofile=None, seqfile=None, seed=seed + 1)
    e3 = pyvolve.Evolver(partitions=p3, tree=pyvolve.read_tree(tree=tree_flank_str))
    e3(ratefile=None, infofile=None, seqfile=None, seed=seed + 2)

    s1, s2, s3 = e1.get_sequences(), e2.get_sequences(), e3.get_sequences()
    all_taxa = sorted(list(s1.keys()))

    if is_complete_ghost:
        # Drop all Clade B taxa from the alignment!
        sampled_taxa = [t for t in all_taxa if not t.startswith("B")]
    else:
        # Drop B0 (the exact donor), leaving only divergent relatives B1, B2, B3
        sampled_taxa = [t for t in all_taxa if t != "B0"]

    seqs = {t: s1[t] + s2[t] + s3[t] for t in sampled_taxa}
    mat = seqs_to_matrix(sampled_taxa, seqs, L)
    meta = {
        "donor_divergence": donor_divergence,
        "is_complete_ghost": is_complete_ghost,
        "recombinant": "Rec",
        "sampled_taxa": sampled_taxa,
        "true_bps": [tract_start, tract_start + tract_len],
        "tract_len": tract_len
    }
    return mat, sampled_taxa, L, meta


# -----------------------------------------------------------------------------
# Option C: Circulating Recombinant Forms (CRFs / Clade Expansion)
# -----------------------------------------------------------------------------

def generate_crf_alignment(
    M_recombinants: int = 4,
    N_total: int = 25,
    divergence: float = 0.05,
    clade_drift: float = 0.005,
    tract_len: int = 600,
    tract_start: int = 1000,
    L: int = 3000,
    seed: int = 42
) -> Tuple[np.ndarray, List[str], int, Dict]:
    """
    An ancestral recombination event occurred at [tract_start, tract_start + tract_len].
    The recombinant lineage expanded into an extant sub-clade of M recombinant taxa (Rec_0, ..., Rec_{M-1}).
    The remaining N_total - M taxa belong to parental Clades A and B.
    """
    np.random.seed(seed)
    d = divergence / 2.0
    n_parents = N_total - M_recombinants
    n_a = n_parents // 2
    n_b = n_parents - n_a

    a_taxa = [f"A{i}" for i in range(n_a)]
    b_taxa = [f"B{i}" for i in range(n_b)]
    rec_taxa = [f"Rec_{i}" for i in range(M_recombinants)]

    c1_str = ",".join(f"{t}:0.0100" for t in a_taxa)
    c2_str = ",".join(f"{t}:0.0100" for t in b_taxa)
    rec_str = ",".join(f"{t}:{clade_drift:.5f}" for t in rec_taxa)

    # Flanks: The CRF clade is sister to Clade A
    tree_flank_str = f"((({rec_str}):0.0100, ({c1_str}):0.0100):{d:.5f}, ({c2_str}):{d:.5f});"
    # Cassette: The CRF clade is sister to Clade B
    tree_cass_str = f"(({c1_str}):{d:.5f}, (({rec_str}):0.0100, ({c2_str}):0.0100):{d:.5f});"

    model = pyvolve.Model("nucleotide", {"parameters": {"kappa": 2.5}})
    p1 = pyvolve.Partition(models=model, size=tract_start)
    p2 = pyvolve.Partition(models=model, size=tract_len)
    p3 = pyvolve.Partition(models=model, size=L - tract_start - tract_len)

    e1 = pyvolve.Evolver(partitions=p1, tree=pyvolve.read_tree(tree=tree_flank_str))
    e1(ratefile=None, infofile=None, seqfile=None, seed=seed)
    e2 = pyvolve.Evolver(partitions=p2, tree=pyvolve.read_tree(tree=tree_cass_str))
    e2(ratefile=None, infofile=None, seqfile=None, seed=seed + 1)
    e3 = pyvolve.Evolver(partitions=p3, tree=pyvolve.read_tree(tree=tree_flank_str))
    e3(ratefile=None, infofile=None, seqfile=None, seed=seed + 2)

    s1, s2, s3 = e1.get_sequences(), e2.get_sequences(), e3.get_sequences()
    taxa = sorted(list(s1.keys()))
    seqs = {t: s1[t] + s2[t] + s3[t] for t in taxa}
    mat = seqs_to_matrix(taxa, seqs, L)
    meta = {
        "M_recombinants": M_recombinants,
        "N_total": N_total,
        "rec_taxa": rec_taxa,
        "true_bps": [tract_start, tract_start + tract_len],
        "tract_len": tract_len,
        "clade_drift": clade_drift
    }
    return mat, taxa, L, meta
