# Dataset Provenance: Olabode et al. (2022) HIV-1 Simulation Benchmark

## Canonical Citation
> **Olabode, A. S., Ng, G. T., Wade, K. E., Salnikov, M., Grant, H. E., Dick, D. W., & Poon, A. F. Y. (2022).**  
> *Revisiting the recombinant history of HIV-1 group M with dynamic network community detection.*  
> **Proceedings of the National Academy of Sciences (PNAS)**, 119(19), e2108815119.  
> DOI: [10.1073/pnas.2108815119](https://doi.org/10.1073/pnas.2108815119)

---

## Dataset Architecture & Experimental Cohorts

This benchmark evaluates full-length continuous-time HIV-1 Group M viral genomes ($L = 9{,}000$ nt, 3,000 codons) across 360 alignments:

### 1. Experiment 1: Post-hoc Reference Benchmark (COMET Design)
- **Taxa:** $N = 16$ reference genomes (4 sequences each from HIV-1/M Subtypes A, B, C, D).
- **Sequence Length:** $L = 9,000\,\text{nt}$.
- **Recombination Scenarios:**
  - 1 Breakpoint: nt 4,500 (Codon 1,500), 10 replicates.
  - 2 Breakpoints: nt 3,000 and 6,000 (Codons 1,000 and 2,000), 10 replicates.
  - 3 Breakpoints: nt 2,250, 4,500, and 6,750 (Codons 750, 1,500, 2,250), 10 replicates.
- **Total Datasets:** 30 alignments.

### 2. Experiment 2: Continuous Time-Scaled Simulations with Gamma Rate Variation
- **Taxa:** $N = 37$ taxa sampled from HIV-1 Group M tree topology (subtypes A--D, F--H, J, K).
- **Evolutionary Parameters:** Total tree length scaled to 3.22 substitutions per codon (643.5 years of divergence under HIV-1 clock of $1.67 \times 10^{-3}$ substitutions/site/year); $\kappa = 8.0$; Gamma site rate variation ($\alpha = 1.5, \beta = 3.0$).
- **Recombination Scenarios:**
  - 1 Breakpoint: 100 independent replicates with breakpoint drawn uniformly across $(0, 9000)$ nt.
  - 2 Breakpoints: 100 independent replicates with 2 uniform breakpoints.
  - 3 Breakpoints: 100 independent replicates with 3 uniform breakpoints.
- **Total Datasets:** 300 alignments.

### 3. Experiment 3: Computational Scaling Cohorts
- **Taxa Count:** $N \in \{50, 100, 200\}$ taxa.
- **Sequence Length:** $L = 9,000\,\text{nt}$.
- **Replicates:** 10 replicates per sample size.
- **Total Datasets:** 30 alignments.

---

## Benchmarked Methods & Baseline Sources
1. **RhizAeon:** Evaluated locally using the calibrated standard workflow:
   `PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=True)`
   `run_recursive_partition_fda_screen(engine, taxa_names=taxa, min_z=1.8, min_pir=0.06, crossover_validation=True, crossover_p_threshold=0.005, min_informative_sites=3, polish_ml=True)`
2. **3SEQ:** Evaluated locally using the official 3SEQ binary (`-full -p 3seq_ptable_250 -q -L40`).
3. **DSBM (Dynamic Stochastic Block Model):** Published distributions from Olabode et al. (2022) Fig. 1A, 1B, and Table S1.
4. **GARD:** Published distributions from Olabode et al. (2022) Fig. 1A, 1B.
5. **RDP4 & RDP5:** Published distributions from Olabode et al. (2022) Fig. 1A, 1B, and Table S1.
