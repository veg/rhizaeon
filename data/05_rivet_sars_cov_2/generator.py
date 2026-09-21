"""
simulations/02_rivet_sars_cov_2/generator.py
============================================
Faithful implementation of the Thornlow & Turakhia SARS-CoV-2 simulation protocol
(bpt26/recombination/makeRandomRecombinants.py; Nature 2022, Bioinformatics 2023).

Simulates recombinant and non-recombinant SARS-CoV-2 whole genomes across the canonical
29,903 nt Wuhan-Hu-1 reference sequence under varying parental divergence (d),
post-recombination mutations (m), and breakpoint topologies (1-BP, 2-BP, Negative Control).
"""

import os
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import numpy as np

DATA_DIR = Path(__file__).resolve().parent / "data"
WUHAN_REF_PATH = DATA_DIR / "wuhan.ref.fa"

NT_MAP = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
NT_CHARS = np.array(["A", "C", "G", "T", "N"], dtype="<U1")


def load_wuhan_reference(ref_path: Optional[Path] = None) -> np.ndarray:
    """Loads and integer-encodes the Wuhan-Hu-1 reference genome (29,903 nt)."""
    p = ref_path or WUHAN_REF_PATH
    if not p.exists():
        raise FileNotFoundError(f"Reference genome not found at {p}")
    seq_chars = []
    with open(p, "r") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq_chars.append(line.upper())
    ref_seq = "".join(seq_chars)
    encoded = np.array([NT_MAP.get(c, 0) for c in ref_seq], dtype=np.int8)
    return encoded


_CACHED_REF: Optional[np.ndarray] = None


def get_wuhan_ref() -> np.ndarray:
    global _CACHED_REF
    if _CACHED_REF is None:
        _CACHED_REF = load_wuhan_reference()
    return _CACHED_REF.copy()


def generate_sars_cov_2_trio(
    n_bps: int,
    d: int,
    m: int,
    seed: int,
    min_segment_len: int = 1500
) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
    """
    Generates a 3-sequence alignment (Recombinant, Donor, Acceptor) of length 29,903 nt.

    Parameters:
      n_bps: 0 (negative control non-recombinant), 1 (single crossover), or 2 (double crossover)
      d: mutational distance (number of substitutions) separating Donor and Acceptor
      m: private post-recombination mutations added to Recombinant
      seed: PRNG seed for exact reproducibility
      min_segment_len: minimum distance from chromosome terminus or between breakpoints

    Returns:
      seq_mat: np.ndarray of shape (3, 29903) with values in {0, 1, 2, 3}
      taxa: ["Recombinant", "Donor", "Acceptor"]
      meta: dictionary with true breakpoints, mutation counts, and coordinates
    """
    rng = np.random.default_rng(seed)
    ref = get_wuhan_ref()
    L = len(ref)

    # 1. Acceptor is derived from reference with background mutations
    acceptor = ref.copy()
    donor = ref.copy()

    # 2. Pick d segregating sites uniformly across the genome to separate Donor and Acceptor
    candidate_sites = np.arange(100, L - 100) # Exclude extreme terminal unsequenced ends
    seg_sites = np.sort(rng.choice(candidate_sites, size=d, replace=False))

    # Mutate Donor at these sites to an alternate nucleotide
    for site in seg_sites:
        orig = donor[site]
        alts = [nt for nt in (0, 1, 2, 3) if nt != orig]
        donor[site] = rng.choice(alts)

    true_bps: List[int] = []
    recombinant = np.zeros(L, dtype=np.int8)

    if n_bps == 0:
        # Negative control: Recombinant is a pure vertical descendant of Acceptor
        recombinant = acceptor.copy()
    elif n_bps == 1:
        # 1-Breakpoint mosaic: Recombinant = Acceptor [:bp] + Donor [bp:]
        bp1 = int(rng.integers(min_segment_len, L - min_segment_len))
        true_bps.append(bp1)
        recombinant[:bp1] = acceptor[:bp1]
        recombinant[bp1:] = donor[bp1:]
    elif n_bps == 2:
        # 2-Breakpoint mosaic: Recombinant = Acceptor [:bp1] + Donor [bp1:bp2] + Acceptor [bp2:]
        bp1 = int(rng.integers(min_segment_len, L - 2 * min_segment_len))
        bp2 = int(rng.integers(bp1 + min_segment_len, L - min_segment_len))
        true_bps = [bp1, bp2]
        recombinant[:bp1] = acceptor[:bp1]
        recombinant[bp1:bp2] = donor[bp1:bp2]
        recombinant[bp2:] = acceptor[bp2:]
    else:
        raise ValueError(f"Unsupported n_bps={n_bps}. Must be 0, 1, or 2.")

    # 3. Add m private post-recombination mutations to Recombinant
    if m > 0:
        available_sites = np.setdiff1d(candidate_sites, seg_sites)
        priv_sites = rng.choice(available_sites, size=m, replace=False)
        for site in priv_sites:
            orig = recombinant[site]
            alts = [nt for nt in (0, 1, 2, 3) if nt != orig]
            recombinant[site] = rng.choice(alts)
    else:
        priv_sites = np.array([], dtype=int)

    seq_mat = np.stack([recombinant, donor, acceptor], axis=0)
    taxa = ["Recombinant", "Donor", "Acceptor"]

    # Compute informative sites supporting the crossover
    # Informative site: Donor and Acceptor differ, Recombinant matches one of them
    diff_mask = (donor != acceptor)
    rec_matches_donor = diff_mask & (recombinant == donor)
    rec_matches_acceptor = diff_mask & (recombinant == acceptor)

    meta = {
        "n_bps": n_bps,
        "d": d,
        "m": m,
        "seed": seed,
        "genome_length": L,
        "true_bps": true_bps,
        "seg_sites_count": int(np.sum(diff_mask)),
        "rec_donor_matches": int(np.sum(rec_matches_donor)),
        "rec_acceptor_matches": int(np.sum(rec_matches_acceptor)),
        "private_mutations_count": int(len(priv_sites))
    }

    return seq_mat, taxa, meta
