"""
benchmarks/didelot_2015_benchmark/generator.py
==============================================
Parametric simulator of bacterial microevolution and homologous gene conversion
replicating the ClonalFrame / ClonalFrameML / SimBac experimental design:
Didelot & Wilson (2015, PLoS Comput Biol 11:e1004041)
Didelot & Falush (2007, Genetics 175:1251-1266)

Key parameters:
- r/m: Ratio of recombination to mutation rates (substitutions from imports vs mutations).
- delta: Mean length of imported gene conversion tracts (nt).
- nu: Mean sequence divergence of imported tracts from donor lineages.
- N: Number of sampled bacterial isolates.
- L: Genome / alignment length (nt).
- theta: Mutation rate per site per unit coalescent time.
"""

import numpy as np
from typing import Dict, List, Tuple, Optional, Any

NT_CHARS = ["A", "C", "G", "T"]
PI = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float64)

def simulate_clonal_tree(n: int = 20, seed: Optional[int] = None) -> Tuple[int, Dict[int, List[int]], Dict[int, float]]:
    """Simulates a Kingman coalescent tree for N bacterial isolates."""
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
        denom = k * (k - 1) / 2.0
        dt = np.random.exponential(1.0 / denom)
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


def simulate_bacterial_gene_conversion(
    n: int = 20,
    L: int = 25000,
    r_over_m: float = 1.0,
    delta: float = 1000.0,
    nu: float = 0.02,
    theta: float = 0.01,
    external_donor_prob: float = 0.5,
    seed: Optional[int] = None
) -> Tuple[List[str], List[Dict[str, Any]], Dict[str, Any]]:
    """
    Simulates bacterial sequence evolution under the ClonalFrame model:
    - Point mutations occur along tree branches at rate theta per site.
    - Recombination events occur along branches at rate R per site, where:
        r/m = (R / theta) * delta * nu  =>  R = (r_over_m * theta) / (delta * nu)
    - Each recombination event imports a tract of mean length delta.
    - If r_over_m == 0.0, the sequence evolves strictly clonally (clonal null control).
    
    Returns:
      taxa_sequences: list of string genomes of length L for N isolates
      imported_tracts: list of metadata dicts for each imported event
      metadata: summary statistics and true breakpoint coordinates
    """
    if seed is not None:
        np.random.seed(seed)
        
    root, children, branch_len = simulate_clonal_tree(n=n, seed=seed)
    
    # Calculate recombination rate per site R
    if r_over_m > 0.0 and delta > 0.0 and nu > 0.0:
        R_rate = (r_over_m * theta) / (delta * nu)
    else:
        R_rate = 0.0
        
    # Sample ancestral root sequence (discrete integers 0, 1, 2, 3)
    root_seq = np.random.choice(4, size=L, p=PI)
    
    node_seqs = {root: root_seq}
    imported_tracts = []
    
    # Top-down tree traversal
    def traverse(node):
        if node not in children:
            return
            
        for ch in children[node]:
            bl = branch_len[ch]
            parent_seq = node_seqs[node]
            curr_seq = np.copy(parent_seq)
            
            # 1. Point mutations: Poisson(theta * bl * L)
            n_muts = np.random.poisson(theta * bl * L)
            if n_muts > 0:
                mut_sites = np.random.choice(L, size=min(n_muts, L), replace=False)
                for s in mut_sites:
                    old_c = curr_seq[s]
                    curr_seq[s] = (old_c + np.random.choice([1, 2, 3])) % 4
                    
            # 2. Homologous gene conversion imports: Poisson(R_rate * bl * L)
            if R_rate > 0.0:
                n_imports = np.random.poisson(R_rate * bl * L)
                for _ in range(n_imports):
                    # Start position
                    start_pos = np.random.randint(0, L - 10)
                    # Tract length ~ Geometric / Exponential(delta)
                    tract_len = max(50, int(np.random.exponential(delta)))
                    end_pos = min(L, start_pos + tract_len)
                    actual_len = end_pos - start_pos
                    
                    # Donor sequence
                    if np.random.rand() < external_donor_prob:
                        donor_type = "external"
                        # Mutate recipient segment at divergence rate nu
                        donor_seg = np.copy(curr_seq[start_pos:end_pos])
                        n_seg_muts = np.random.poisson(nu * actual_len)
                        if n_seg_muts > 0:
                            seg_sites = np.random.choice(actual_len, size=min(n_seg_muts, actual_len), replace=False)
                            for ss in seg_sites:
                                donor_seg[ss] = (donor_seg[ss] + np.random.choice([1, 2, 3])) % 4
                    else:
                        donor_type = "internal"
                        # Borrow from root or another node
                        donor_node = np.random.choice(list(node_seqs.keys()))
                        donor_seg = np.copy(node_seqs[donor_node][start_pos:end_pos])
                        # Apply subtle divergence
                        n_seg_muts = np.random.poisson(nu * 0.5 * actual_len)
                        if n_seg_muts > 0:
                            seg_sites = np.random.choice(actual_len, size=min(n_seg_muts, actual_len), replace=False)
                            for ss in seg_sites:
                                donor_seg[ss] = (donor_seg[ss] + np.random.choice([1, 2, 3])) % 4
                                
                    # Replace tract
                    curr_seq[start_pos:end_pos] = donor_seg
                    
                    imported_tracts.append({
                        "recipient_node": ch,
                        "is_leaf": (ch < n),
                        "start": int(start_pos),
                        "end": int(end_pos),
                        "length": int(actual_len),
                        "donor_type": donor_type,
                        "branch_len": float(bl)
                    })
                    
            node_seqs[ch] = curr_seq
            traverse(ch)
            
    traverse(root)
    
    # Extract leaf sequences
    taxa_seqs = ["".join(NT_CHARS[c] for c in node_seqs[i]) for i in range(n)]
    
    # True breakpoints are the boundaries of imported tracts
    all_bps = set()
    for tr in imported_tracts:
        if tr["start"] > 0:
            all_bps.add(tr["start"])
        if tr["end"] < L:
            all_bps.add(tr["end"])
            
    sorted_bps = sorted(all_bps)
    
    meta = {
        "n_taxa": n,
        "length": L,
        "r_over_m": r_over_m,
        "delta": delta,
        "nu": nu,
        "theta": theta,
        "num_imports": len(imported_tracts),
        "num_true_bps": len(sorted_bps),
        "true_breakpoints": sorted_bps,
        "is_recombinant": (len(imported_tracts) > 0 and len(sorted_bps) > 0)
    }
    
    return taxa_seqs, imported_tracts, meta
