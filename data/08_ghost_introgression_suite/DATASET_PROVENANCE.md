# Dataset Provenance: Ghost Introgression & Deep Phylogenies Benchmark Suite

## Canonical Citations & Conceptual Framework
> **Kosakovsky Pond, S. L., et al. (2026).**  
> *Clean Breaks: Tree-Free Discovery of Mosaic Genomes on Continuous Sequence Manifolds.*  
> In preparation.

> **Martin, D. P., et al. (2021).**  
> *RDP5: A computer program for analysing recombination in, and removing recombination from, multiple sequence alignments.*  
> **Molecular Biology and Evolution**, 38(6), 2690–2697.  
> DOI: [10.1093/molbev/msab063](https://doi.org/10.1093/molbev/msab063)

> **Boni, M. F., Posada, D., & Feldman, M. W. (2007).**  
> *An exact nonparametric method for inferring mosaic structure with application to genetics data and simulated data.*  
> **Genetics**, 176(2), 1035–1047.  
> DOI: [10.1534/genetics.106.068874](https://doi.org/10.1534/genetics.106.068874)

---

## Benchmark Design & Experimental Conditions

Natural pathogen transmission and evolutionary history frequently lack direct parental donors due to incomplete surveillance, extinct reservoir lineages, or deep evolutionary divergence. This benchmark evaluates recombination detection across 470 synthetic multi-taxon alignments ($L = 3{,}000$ nt) spanning two primary axes:

### Axis 1: Ghost Lineage Introgression and Donor Drift ($N = 25$ Taxa, 120 Alignments)
- Recombinant sequence $R$ inherits an ancestral 600-nt cassette from a donor lineage that underwent independent post-recombination divergence ($\Delta d$) before sampling:
  1. $\Delta d = 0.00$ (Exact sampled donor, baseline)
  2. $\Delta d = 0.01$ (Minimal drift, 1% divergence)
  3. $\Delta d = 0.03$ (Moderate drift, 3% divergence)
  4. $\Delta d = 0.05$ (Substantial drift, 5% divergence)
  5. $\Delta d = 0.10$ (Deep divergence, 10% drift)
  6. 100% Clade Missing (Complete ghost introgression: all taxa belonging to the donor clade are excised prior to analysis).
- 20 independent replicates per parameter condition (120 alignments total).

### Axis 2: Deep Taxonomic Scaling & Clonal Null Calibration ($N \in \{5, 10, 20, 50, 100\}$, 250 Alignments)
- Evaluates detection power, dual-breakpoint resolution, and false positive rate under expanding phylogenetic cohorts ($N = 5, 10, 20, 50, 100$) at low ($d = 0.02$) and moderate ($d = 0.05$) sequence divergence.
- 10 non-recombinant clonal null replicates per condition (evaluating Type I error control against parental shadowing).
- 15 recombinant replicates per condition (3,000 nt, cassette boundaries at nt 1,000 and nt 1,500).

---

## Standardized RhizAeon Execution Workflow
Evaluations executed under the calibrated standard RhizAeon pipeline:
```python
engine = PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=True)
bps = run_recursive_partition_fda_screen(
    engine,
    taxa_names=taxa_names,
    min_z=1.8,
    min_pir=0.06,
    crossover_validation=True,
    crossover_p_threshold=0.005,
    min_informative_sites=3,
    polish_ml=True
)
```
