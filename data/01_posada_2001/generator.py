"""
benchmarks/posada_2001_benchmark/generator.py
=============================================
Parametric coalescent simulator replicating the exact experimental design of
Posada & Crandall (2001, PNAS 98:13757-13762):
- Coalescent with recombination (Hudson 1983; Kaplan & Hudson 1985)
- Sample size n = 10 taxa, sequence length L = 1,000 nt
- Mutation parameter theta in {10, 50, 100, 200}
- Recombination parameter rho in {0, 1, 4, 16, 64}
- Gamma rate variation alpha in {None (inf), 2.0, 0.5, 0.05}
- Nucleotide substitution model: HKY85 (kappa = 2.0, pi = [0.4, 0.2, 0.3, 0.1])
"""

import numpy as np
from scipy.linalg import expm
from typing import Dict, List, Tuple, Optional, Any

NT_CHARS = ["A", "C", "G", "T"]

# HKY85 equilibrium frequencies: A=0.40, C=0.20, G=0.30, T=0.10 (Posada & Crandall 2001 Table 4)
# Note: In Table 4, order is pi_A=0.40, pi_C=0.20, pi_T=0.10, pi_G=0.30
PI = np.array([0.40, 0.20, 0.30, 0.10], dtype=np.float64)
KAPPA = 2.0

def build_hky_rate_matrix(pi: np.ndarray = PI, kappa: float = KAPPA) -> np.ndarray:
    """Constructs normalized continuous-time HKY generator Q matrix."""
    Q = np.zeros((4, 4), dtype=np.float64)
    # State mapping: 0: A, 1: C, 2: G, 3: T
    # Purines: 0 (A), 2 (G); Pyrimidines: 1 (C), 3 (T)
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            is_transition = (i in (0, 2) and j in (0, 2)) or (i in (1, 3) and j in (1, 3))
            mult = kappa if is_transition else 1.0
            Q[i, j] = mult * pi[j]
            
    for i in range(4):
        Q[i, i] = -np.sum(Q[i, :])
        
    # Scale Q so that average substitution rate -sum(pi_i * Q_ii) = 1.0
    average_rate = -np.sum(pi * np.diag(Q))
    Q /= average_rate
    return Q

GLOBAL_HKY_Q = build_hky_rate_matrix()


def simulate_pure_coalescent_tree(n: int = 10, seed: Optional[int] = None) -> Tuple[int, Dict[int, List[int]], Dict[int, float]]:
    """
    Simulates a standard neutral Kingman coalescent tree with n leaves (rho = 0).
    Returns (root_id, children_dict, branch_length_dict).
    Branch lengths are in units of 2N generations.
    """
    if seed is not None:
        np.random.seed(seed)
        
    active = list(range(n))
    node_time = {i: 0.0 for i in range(n)}
    children = {}
    branch_len = {}
    next_node = n
    t = 0.0
    
    while len(active) > 1:
        k = len(active)
        rate = k * (k - 1) / 2.0
        dt = np.random.exponential(1.0 / rate)
        t += dt
        
        pair_idx = np.random.choice(len(active), size=2, replace=False)
        u, v = active[pair_idx[0]], active[pair_idx[1]]
        
        node_time[next_node] = t
        children[next_node] = [u, v]
        branch_len[u] = t - node_time[u]
        branch_len[v] = t - node_time[v]
        
        active.remove(u)
        active.remove(v)
        active.append(next_node)
        next_node += 1
        
    root = active[0]
    return root, children, branch_len


def simulate_hudson_arg(
    n: int = 10,
    L: int = 1000,
    rho: float = 4.0,
    seed: Optional[int] = None
) -> List[Tuple[int, int, int, Dict[int, List[int]], Dict[int, float]]]:
    """
    Simulates the ancestral recombination graph (ARG) backwards in time following
    Hudson (1983) and Kaplan & Hudson (1985) using high-performance vectorized site masks.
    """
    if rho <= 0.0:
        root, ch, bl = simulate_pure_coalescent_tree(n=n, seed=seed)
        return [(0, L, root, ch, bl)]
        
    if seed is not None:
        np.random.seed(seed)
        
    lineages = {i: np.ones(L, dtype=bool) for i in range(n)}
    node_time = {i: 0.0 for i in range(n)}
    active_node = {i: i for i in range(n)}
    site_parents = {s: {} for s in range(L)}
    next_node = n
    
    site_k = np.full(L, n, dtype=np.int32)
    t = 0.0
    
    max_steps = 50000
    step = 0
    
    while np.any(site_k > 1) and step < max_steps:
        step += 1
        g_vals = {}
        for lid, mask in lineages.items():
            true_idx = np.where(mask)[0]
            if len(true_idx) > 1:
                g_vals[lid] = int(true_idx[-1] - true_idx[0])
            else:
                g_vals[lid] = 0
                
        G = sum(g_vals.values())
        k = len(lineages)
        if k <= 1 and G == 0:
            break
            
        lambda_coal = k * (k - 1) / 2.0
        lambda_rec = (rho * G) / (2.0 * L) if (L > 0 and rho > 0 and G > 0) else 0.0
        tot = lambda_coal + lambda_rec
        
        dt = np.random.exponential(1.0 / tot)
        t += dt
        
        if np.random.rand() < (lambda_coal / tot):
            # Coalescence
            lids = list(lineages.keys())
            idx = np.random.choice(len(lids), size=2, replace=False)
            u, v = lids[idx[0]], lids[idx[1]]
            m_u, m_v = lineages[u], lineages[v]
            overlap = m_u & m_v
            
            new_node = next_node
            next_node += 1
            node_time[new_node] = t
            
            node_u = active_node[u]
            node_v = active_node[v]
            
            for s in np.where(overlap)[0]:
                site_parents[s][node_u] = (new_node, t - node_time[node_u])
                site_parents[s][node_v] = (new_node, t - node_time[node_v])
            site_k[overlap] -= 1
            
            for s in np.where(m_u & ~overlap)[0]:
                site_parents[s][node_u] = (new_node, t - node_time[node_u])
            for s in np.where(m_v & ~overlap)[0]:
                site_parents[s][node_v] = (new_node, t - node_time[node_v])
                
            del lineages[u]
            del lineages[v]
            del active_node[u]
            del active_node[v]
            
            union_active = (m_u | m_v) & (site_k > 1)
            if np.any(union_active):
                lineages[new_node] = union_active
                active_node[new_node] = new_node
        else:
            # Recombination
            eligible = [lid for lid in lineages if g_vals[lid] > 0]
            w = [g_vals[lid] for lid in eligible]
            lid = eligible[np.random.choice(len(eligible), p=np.array(w) / sum(w))]
            mask = lineages[lid]
            true_idx = np.where(mask)[0]
            b = np.random.randint(true_idx[0], true_idx[-1] + 1)
            
            m_l = mask.copy()
            m_l[b + 1:] = False
            m_r = mask.copy()
            m_r[:b + 1] = False
            
            old_node = active_node[lid]
            del lineages[lid]
            del active_node[lid]
            
            if np.any(m_l):
                nl = next_node
                next_node += 1
                node_time[nl] = t
                lineages[nl] = m_l
                active_node[nl] = nl
                for s in np.where(m_l)[0]:
                    site_parents[s][old_node] = (nl, t - node_time[old_node])
            if np.any(m_r):
                nr = next_node
                next_node += 1
                node_time[nr] = t
                lineages[nr] = m_r
                active_node[nr] = nr
                for s in np.where(m_r)[0]:
                    site_parents[s][old_node] = (nr, t - node_time[old_node])

    # Convert per-site parent maps into contiguous partitions
    # To compare site trees, we hash the sorted child -> (parent, rounded_blen) mapping
    def tree_key(s):
        return tuple(sorted((c, p[0], round(p[1], 6)) for c, p in site_parents[s].items()))
        
    partitions = []
    curr_start = 0
    curr_key = tree_key(0)
    for s in range(1, L):
        k = tree_key(s)
        if k != curr_key:
            partitions.append((curr_start, s))
            curr_start = s
            curr_key = k
    partitions.append((curr_start, L))
    
    # For each partition, construct the rooted tree graph (root, children, branch_len)
    result = []
    for start, end in partitions:
        mid_site = start
        sp = site_parents[mid_site]
        # Build tree representation for leaves 0..n-1
        children = {}
        branch_len = {}
        # Find root: the node that is a parent but has no parent itself
        all_children = set(sp.keys())
        all_parents = set(p[0] for p in sp.values())
        root_candidates = all_parents - all_children
        root = list(root_candidates)[0] if root_candidates else max(all_parents)
        
        for c, (p, bl) in sp.items():
            children.setdefault(p, []).append(c)
            branch_len[c] = bl
            
        result.append((start, end, root, children, branch_len))
        
    return result


def evolve_sequences(
    partitions: List[Tuple[int, int, int, Dict[int, List[int]], Dict[int, float]]],
    n: int = 10,
    L: int = 1000,
    theta: float = 50.0,
    alpha: Optional[float] = None,
    seed: Optional[int] = None
) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
    """
    Evolves nucleotide sequences along the partitioned trees under HKY85 with optional
    Gamma rate heterogeneity across sites.
    
    Returns:
      (matrix [n, L], taxa_names, metadata)
    """
    if seed is not None:
        np.random.seed(seed)
        
    # Scale factor from 2N coalescent time to expected substitutions per site
    # theta = 4N * mu * L => 2N * mu = theta / (2 * L)
    mu_base = theta / (2.0 * L)
    
    # Among-site rate variation: r_s ~ Gamma(alpha, 1/alpha) so E[r_s] = 1.0
    if alpha is not None and alpha > 0:
        site_rates = np.random.gamma(shape=alpha, scale=1.0 / alpha, size=L)
    else:
        site_rates = np.ones(L, dtype=np.float64)
        
    taxa = [f"Taxon_{i}" for i in range(n)]
    mat = np.zeros((n, L), dtype=np.int8)
    
    for start, end, root, children, branch_len in partitions:
        seg_len = end - start
        if seg_len <= 0:
            continue
            
        seg_rates = site_rates[start:end]
        
        # Sample root state
        root_seq = np.random.choice(4, size=seg_len, p=PI)
        node_seqs = {root: root_seq}
        
        # Preorder traversal from root to leaves
        stack = [root]
        while stack:
            curr = stack.pop()
            for ch in children.get(curr, []):
                bl_raw = branch_len.get(ch, 0.0)
                curr_seq = node_seqs[curr]
                ch_seq = np.zeros(seg_len, dtype=np.int8)
                
                # If no rate variation, P matrix is constant for this branch
                if alpha is None:
                    eff_bl = bl_raw * mu_base
                    if eff_bl > 50.0:
                        P = np.tile(PI, (4, 1))
                    else:
                        P = expm(GLOBAL_HKY_Q * eff_bl)
                        P = np.nan_to_num(P, nan=0.25)
                        P = np.clip(P, 0.0, 1.0)
                        row_sums = np.sum(P, axis=1, keepdims=True)
                        row_sums[row_sums == 0] = 1.0
                        P /= row_sums
                    for s in range(4):
                        mask = (curr_seq == s)
                        if np.any(mask):
                            ch_seq[mask] = np.random.choice(4, size=np.sum(mask), p=P[s])
                else:
                    # Site-specific rates: group sites by approximate rate or loop
                    eff_bls = bl_raw * mu_base * seg_rates
                    unique_rates, inv = np.unique(np.round(eff_bls, 4), return_inverse=True)
                    for idx_u, r in enumerate(unique_rates):
                        site_mask = (inv == idx_u)
                        if not np.any(site_mask):
                            continue
                        if r > 50.0:
                            P = np.tile(PI, (4, 1))
                        else:
                            P = expm(GLOBAL_HKY_Q * r)
                            P = np.nan_to_num(P, nan=0.25)
                            P = np.clip(P, 0.0, 1.0)
                            row_sums = np.sum(P, axis=1, keepdims=True)
                            row_sums[row_sums == 0] = 1.0
                            P /= row_sums
                        sub_curr = curr_seq[site_mask]
                        sub_ch = np.zeros(len(sub_curr), dtype=np.int8)
                        for s in range(4):
                            m = (sub_curr == s)
                            if np.any(m):
                                sub_ch[m] = np.random.choice(4, size=np.sum(m), p=P[s])
                        ch_seq[site_mask] = sub_ch
                        
                node_seqs[ch] = ch_seq
                stack.append(ch)
                
        # Fill in leaves 0..n-1
        for i in range(n):
            mat[i, start:end] = node_seqs[i]
            
    # Calculate true breakpoints
    true_breakpoints = [p[1] for p in partitions[:-1]]
    
    metadata = {
        "n": n,
        "L": L,
        "theta": theta,
        "alpha": alpha,
        "num_partitions": len(partitions),
        "true_breakpoints": true_breakpoints,
        "is_recombinant": len(partitions) > 1
    }
    
    return mat, taxa, metadata


def generate_posada_dataset(
    n: int = 10,
    L: int = 1000,
    theta: float = 50.0,
    rho: float = 4.0,
    alpha: Optional[float] = None,
    seed: Optional[int] = None
) -> Tuple[np.ndarray, List[str], Dict[str, Any]]:
    """
    Master generator function for Posada & Crandall (2001) replication:
    Simulates ARG and evolves HKY85 sequences under specified parameters.
    """
    partitions = simulate_hudson_arg(n=n, L=L, rho=rho, seed=seed)
    # Use derived seed for evolution
    evo_seed = (seed + 1000003) if seed is not None else None
    mat, taxa, meta = evolve_sequences(
        partitions=partitions,
        n=n,
        L=L,
        theta=theta,
        alpha=alpha,
        seed=evo_seed
    )
    meta["rho"] = rho
    return mat, taxa, meta
