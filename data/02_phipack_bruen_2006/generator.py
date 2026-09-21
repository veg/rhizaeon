"""
simulations/02_phipack_bruen_2006/generator.py
==============================================
Parametric coalescent simulator replicating the experimental framework of
Bruen, Philippe, and Bryant (2006, Genetics 172:2665-2681):
- Exact Ancestral Recombination Graph (Hudson 1983; Kaplan & Hudson 1985 via msprime)
- Variable sample size n in {10, 50} taxa; sequence length L = 1,000 nt
- Population mutation parameter theta = 4N*mu in {5, 10, 20}
- Recombination intensity rho = 4N*r in {0, 1, 4, 16, 64}
- Demographic growth: Exponential population growth (rate beta in {0, 5, 20}, star-like genealogies)
- Site-to-site rate variation: Continuous Gamma distribution with shape parameter alpha in {0.5, 2.0}
- Continuous-time Markov substitution under HKY85 (kappa = 2.0, equal base frequencies)
"""

import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"

import msprime
import numpy as np
from typing import Dict, List, Tuple, Optional, Any

NT_CHARS = ["A", "C", "G", "T"]
BASE_FREQS = [0.25, 0.25, 0.25, 0.25]
DEFAULT_KAPPA = 2.0
REFERENCE_N0 = 10000.0


def generate_phipack_dataset(
    n: int = 10,
    L: int = 1000,
    theta: float = 10.0,
    rho: float = 0.0,
    beta: float = 0.0,
    alpha: Optional[float] = None,
    seed: Optional[int] = None,
    kappa: float = DEFAULT_KAPPA
) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
    """
    Simulates a synthetic alignment under the Bruen et al. (2006) protocol.
    
    Parameters
    ----------
    n : int
        Number of taxa (monoploid sample sequences). Typically 10 or 50.
    L : int
        Sequence length in nucleotides (1,000 nt).
    theta : float
        Population mutation rate theta = 4N0 * mu_locus.
    rho : float
        Population recombination rate rho = 4N0 * r_locus.
    beta : float
        Exponential growth rate parameter (0.0 for constant size, 5.0 or 20.0 for expansion).
        Backwards in time, N(t) = N0 * exp(-g * t) where g = beta / (4 * N0).
    alpha : float or None
        Gamma shape parameter for site rate variation. None implies uniform rates.
    seed : int or None
        Random seed for reproducibility.
    kappa : float
        Transition/transversion ratio for HKY85 (default 2.0).

    Returns
    -------
    msa : np.ndarray of shape (n, L) with values in {0, 1, 2, 3}
    taxa : List[str]
    meta : Dict[str, Any]
    """
    N0 = REFERENCE_N0
    g = (beta / (4.0 * N0)) if beta > 0.0 else 0.0

    demography = None
    if g > 0.0:
        demography = msprime.Demography()
        demography.add_population(name="pop", initial_size=N0, growth_rate=g)

    # Convert total locus parameters to per-base rates in generations
    recomb_rate = (rho / (4.0 * N0 * float(max(1, L - 1)))) if rho > 0.0 else 0.0
    base_mu = theta / (4.0 * N0 * float(L))

    # 1. Simulate Ancestral Recombination Graph (ARG)
    if demography is not None:
        ts = msprime.sim_ancestry(
            samples=[msprime.SampleSet(n, ploidy=1)],
            sequence_length=L,
            recombination_rate=recomb_rate,
            demography=demography,
            random_seed=seed,
            record_provenance=False
        )
    else:
        ts = msprime.sim_ancestry(
            samples=[msprime.SampleSet(n, ploidy=1)],
            sequence_length=L,
            recombination_rate=recomb_rate,
            population_size=N0,
            random_seed=seed,
            record_provenance=False
        )

    # Extract true breakpoint positions from marginal tree boundaries
    true_bps = [int(round(t.interval.right)) for t in ts.trees()][:-1]

    # 2. Setup site-to-site mutation rate variation (Gamma or uniform)
    rng = np.random.default_rng(seed)
    if alpha is not None and alpha > 0.0:
        gamma_multipliers = rng.gamma(shape=alpha, scale=1.0 / alpha, size=L)
        gamma_multipliers /= np.mean(gamma_multipliers)
        rate_map = msprime.RateMap(
            position=list(range(L + 1)),
            rate=list(gamma_multipliers * base_mu)
        )
        effective_rate = rate_map
    else:
        effective_rate = base_mu

    # 3. Simulate finite-sites substitutions along ARG under HKY85
    hky_model = msprime.HKY(kappa=kappa)
    mts = msprime.sim_mutations(
        ts,
        rate=effective_rate,
        model=hky_model,
        random_seed=seed,
        record_provenance=False
    )

    # 4. Reconstruct full nucleotide alignment matrix (n x L)
    root_seq = rng.choice(4, size=L, p=BASE_FREQS)
    msa = np.tile(root_seq, (n, 1)).astype(np.int8)

    for var in mts.variants():
        pos = int(round(var.site.position))
        if 0 <= pos < L:
            for s_idx, g_idx in enumerate(var.genotypes):
                ch = var.alleles[g_idx]
                if ch in NT_CHARS:
                    msa[s_idx, pos] = NT_CHARS.index(ch)

    taxa = [f"seq_{i}" for i in range(n)]

    # Compute descriptive summary metrics
    diffs = [len(np.unique(msa[:, j])) for j in range(L)]
    seg_sites = int(sum(d > 1 for d in diffs))
    inf_sites = 0
    for j in range(L):
        col = msa[:, j]
        cnt = np.bincount(col, minlength=4)
        if np.sum(cnt >= 2) >= 2:
            inf_sites += 1

    meta = {
        "n": n,
        "L": L,
        "theta": theta,
        "rho": rho,
        "beta": beta,
        "alpha": alpha if alpha is not None else "inf",
        "num_trees": ts.num_trees,
        "true_bps": true_bps,
        "num_bps": len(true_bps),
        "segregating_sites": seg_sites,
        "informative_sites": inf_sites,
        "is_recombinant": (rho > 0.0 and len(true_bps) > 0)
    }

    return msa, taxa, meta
