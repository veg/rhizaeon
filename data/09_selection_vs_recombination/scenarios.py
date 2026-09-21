"""
benchmarks/grand_1000_benchmark/scenarios.py
============================================
Comprehensive 20-Scenario Experimental Design for the Grand-Scale 1,000-Alignment
Recombination Benchmark Suite.

20 Scenarios x 50 Stochastic Replicates = 1,000 Datasets:
  Dimension A: Null Controls & Rate Heterogeneity (FPR & Confounder Resistance)
    1. null_homogeneous: Clock-like homogeneous codon substitution (alpha=inf)
    2. null_srv_mild: Moderate spatial rate variation (alpha=1.0)
    3. null_srv_strong: Severe spatial rate variation (alpha=0.3)
    4. null_selection: Purifying/positive selection heterogeneity (omega=0.85)

  Dimension B: Divergence & Mutation Density Spectrum (Detection Thresholds)
    5. div_ultra_shallow: Scale = 0.05 (~99.5% seq identity; acute outbreak regime)
    6. div_shallow: Scale = 0.20 (~98% identity)
    7. div_medium: Scale = 1.00 (~85% identity; standard viral family)
    8. div_deep: Scale = 2.50 (~65% identity; high divergence)
    9. div_ultra_deep: Scale = 4.00 (~50% identity; saturated divergence)

  Dimension C: Recombination Tract Architecture (Spatial Resolution Limits)
   10. tract_micro_30: Micro-recombinant tract of 30 codons / 90 nt [285, 315]
   11. tract_short_60: Short tract of 60 codons / 180 nt [270, 330]
   12. tract_asymmetric_head: 20% head tract [120] (codons 1-120)
   13. tract_asymmetric_tail: 20% tail tract [480] (codons 480-600)
   14. double_cassette: Cassette insertion / gene conversion [200, 400]
   15. triple_mosaic: Multi-way 3-breakpoint mosaic [150, 300, 450]
   16. quad_mosaic: Multi-way 4-breakpoint mosaic [120, 240, 360, 480]

  Dimension D: Topological Difficulty & Ghost Donors (Phylogenetic Sensitivity)
   17. subtle_intra_clade: Intra-clade SPR move (RF=4, codon 300)
   18. deep_inter_clade: Inter-clade SPR move (RF=11, codon 300)
   19. ghost_donor_close: Unsampled sister clade donor (pruned donor, delta=0.08)
   20. ghost_donor_distant: Unsampled deep lineage donor (pruned donor, delta=0.25)
"""

import re
from typing import Dict, List, Tuple, Any, Optional

# Base Clade Definitions (16 Tips + Ghost nodes)
CLADE_A_BASE = "((A1:0.08,R_mosaic:0.08):0.04,(A2:0.08,A3:0.08):0.04):0.08"
CLADE_A_CLEAN = "(A1:0.08,(A2:0.08,A3:0.08):0.04):0.08"
CLADE_A_SUBTLE = "((A1:0.08,A2:0.08):0.04,(A3:0.08,R_mosaic:0.08):0.04):0.08"

CLADE_B_BASE = "((B1:0.08,B2:0.08):0.04,(B3:0.08,B4:0.08):0.04):0.08"
CLADE_B_RECOMB = "(((B1:0.08,R_mosaic:0.08):0.04,B2:0.08):0.04,(B3:0.08,B4:0.08):0.04):0.08"

CLADE_C_BASE = "((C1:0.08,C2:0.08):0.04,(C3:0.08,C4:0.08):0.04):0.08"
CLADE_C_RECOMB = "(((C1:0.08,R_mosaic:0.08):0.04,C2:0.08):0.04,(C3:0.08,C4:0.08):0.04):0.08"

CLADE_D_BASE = "((D1:0.08,D2:0.08):0.04,(D3:0.08,D4:0.08):0.04):0.08"

# Ghost Clades (with unsampled donor G_donor to be pruned)
CLADE_C_GHOST_CLOSE = "(((C1:0.08,(R_mosaic:0.04,G_donor:0.04):0.04):0.04,C2:0.08):0.04,(C3:0.08,C4:0.08):0.04):0.08"
CLADE_D_GHOST_DISTANT = f"((R_mosaic:0.20,G_donor:0.20):0.10,{CLADE_D_BASE}):0.08"

# 4 Fully Binary Rooted Topologies
T1_BASE = f"(({CLADE_A_BASE},{CLADE_B_BASE}):0.08,({CLADE_C_BASE},{CLADE_D_BASE}):0.08);"
T2_DEEP = f"(({CLADE_A_CLEAN},{CLADE_B_BASE}):0.08,({CLADE_C_RECOMB},{CLADE_D_BASE}):0.08);"
T_SUBTLE = f"(({CLADE_A_SUBTLE},{CLADE_B_BASE}):0.08,({CLADE_C_BASE},{CLADE_D_BASE}):0.08);"
T3_CLADE_B = f"(({CLADE_A_CLEAN},{CLADE_B_RECOMB}):0.08,({CLADE_C_BASE},{CLADE_D_BASE}):0.08);"
T_GHOST_CLOSE = f"(({CLADE_A_CLEAN},{CLADE_B_BASE}):0.08,({CLADE_C_GHOST_CLOSE},{CLADE_D_BASE}):0.08);"
T_GHOST_DIST = f"(({CLADE_A_CLEAN},{CLADE_B_BASE}):0.08,({CLADE_C_BASE},{CLADE_D_GHOST_DISTANT}):0.08);"


def scale_newick(t_str: str, scale: float) -> str:
    """Scales branch lengths by a multiplicative factor."""
    clean = t_str.replace(" ", "").strip()
    return re.sub(r':([0-9\.]+)', lambda m: f':{float(m.group(1)) * scale:.5f}', clean)


GRAND_SCENARIOS: Dict[str, Dict[str, Any]] = {
    # Dimension A: Null Controls & Rate Variation
    "null_homogeneous": {
        "dimension": "Dimension A: Null Controls",
        "description": "Negative control: Clock-like homogeneous substitution rates",
        "partitions": [(600, T1_BASE)],
        "true_bps": [],
        "recombinant_taxa": [],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "null_srv_mild": {
        "dimension": "Dimension A: Null Controls",
        "description": "Negative control: Moderate spatial rate variation (Gamma alpha=1.0)",
        "partitions": [(600, T1_BASE)],
        "true_bps": [],
        "recombinant_taxa": [],
        "srv_alpha": 1.0,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "null_srv_strong": {
        "dimension": "Dimension A: Null Controls",
        "description": "Negative control: Severe spatial rate variation (Gamma alpha=0.3)",
        "partitions": [(600, T1_BASE)],
        "true_bps": [],
        "recombinant_taxa": [],
        "srv_alpha": 0.3,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "null_selection": {
        "dimension": "Dimension A: Null Controls",
        "description": "Negative control: Selection rate variation (omega=0.85)",
        "partitions": [(600, T1_BASE)],
        "true_bps": [],
        "recombinant_taxa": [],
        "srv_alpha": None,
        "omega": 0.85,
        "tree_scale": 1.0,
        "prune_taxa": []
    },

    # Dimension B: Divergence Phase Transitions
    "div_ultra_shallow": {
        "dimension": "Dimension B: Divergence Spectrum",
        "description": "Ultra-low divergence outbreak regime (Scale=0.05, ~99.5% seq identity)",
        "partitions": [(300, T1_BASE), (300, T2_DEEP)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 0.05,
        "prune_taxa": []
    },
    "div_shallow": {
        "dimension": "Dimension B: Divergence Spectrum",
        "description": "Low divergence epidemic regime (Scale=0.20, ~98% seq identity)",
        "partitions": [(300, T1_BASE), (300, T2_DEEP)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 0.20,
        "prune_taxa": []
    },
    "div_medium": {
        "dimension": "Dimension B: Divergence Spectrum",
        "description": "Standard viral divergence regime (Scale=1.00, ~85% seq identity)",
        "partitions": [(300, T1_BASE), (300, T2_DEEP)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.00,
        "prune_taxa": []
    },
    "div_deep": {
        "dimension": "Dimension B: Divergence Spectrum",
        "description": "High divergence inter-clade regime (Scale=2.50, ~65% seq identity)",
        "partitions": [(300, T1_BASE), (300, T2_DEEP)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 2.50,
        "prune_taxa": []
    },
    "div_ultra_deep": {
        "dimension": "Dimension B: Divergence Spectrum",
        "description": "Mutational saturation regime (Scale=4.00, ~50% seq identity)",
        "partitions": [(300, T1_BASE), (300, T2_DEEP)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 4.00,
        "prune_taxa": []
    },

    # Dimension C: Recombination Tract Architecture
    "tract_micro_30": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Micro-recombinant tract (30 codons / 90 nt at [285, 315])",
        "partitions": [(285, T1_BASE), (30, T2_DEEP), (285, T1_BASE)],
        "true_bps": [285, 315],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "tract_short_60": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Short recombinant tract (60 codons / 180 nt at [270, 330])",
        "partitions": [(270, T1_BASE), (60, T2_DEEP), (270, T1_BASE)],
        "true_bps": [270, 330],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "tract_asymmetric_head": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Asymmetric 20% head tract at codon 120",
        "partitions": [(120, T2_DEEP), (480, T1_BASE)],
        "true_bps": [120],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "tract_asymmetric_tail": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Asymmetric 20% tail tract at codon 480",
        "partitions": [(480, T1_BASE), (120, T2_DEEP)],
        "true_bps": [480],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "double_cassette": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Gene conversion cassette insertion at [200, 400]",
        "partitions": [(200, T1_BASE), (200, T2_DEEP), (200, T1_BASE)],
        "true_bps": [200, 400],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "triple_mosaic": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Three-way mosaic reticulation at [150, 300, 450] (T1 -> T2 -> T3 -> T1)",
        "partitions": [(150, T1_BASE), (150, T2_DEEP), (150, T3_CLADE_B), (150, T1_BASE)],
        "true_bps": [150, 300, 450],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "quad_mosaic": {
        "dimension": "Dimension C: Tract Architecture",
        "description": "Four-way mosaic reticulation at [120, 240, 360, 480]",
        "partitions": [(120, T1_BASE), (120, T2_DEEP), (120, T3_CLADE_B), (120, T2_DEEP), (120, T1_BASE)],
        "true_bps": [120, 240, 360, 480],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },

    # Dimension D: Topological Difficulty & Ghost Donors
    "subtle_intra_clade": {
        "dimension": "Dimension D: Topological Difficulty",
        "description": "Subtle intra-clade SPR move (RF=4 at codon 300)",
        "partitions": [(300, T1_BASE), (300, T_SUBTLE)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "deep_inter_clade": {
        "dimension": "Dimension D: Topological Difficulty",
        "description": "Deep inter-clade SPR move (RF=11 at codon 300)",
        "partitions": [(300, T1_BASE), (300, T2_DEEP)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": []
    },
    "ghost_donor_close": {
        "dimension": "Dimension D: Topological Difficulty",
        "description": "Unsampled sister clade ghost donor (G_donor pruned from alignment)",
        "partitions": [(300, T1_BASE), (300, T_GHOST_CLOSE)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 1.0,
        "prune_taxa": ["G_donor"]
    },
    "ghost_donor_distant": {
        "dimension": "Dimension D: Topological Difficulty",
        "description": "Unsampled deep outgroup ghost donor (G_donor pruned, Scale=2.0)",
        "partitions": [(300, T1_BASE), (300, T_GHOST_DIST)],
        "true_bps": [300],
        "recombinant_taxa": ["R_mosaic"],
        "srv_alpha": None,
        "omega": 0.25,
        "tree_scale": 2.0,
        "prune_taxa": ["G_donor"]
    }
}
