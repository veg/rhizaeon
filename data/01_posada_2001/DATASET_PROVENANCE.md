# Dataset Provenance: Posada & Crandall (2001 PNAS)

## Reference Publication
- **Title:** Evaluation of methods for detecting recombination from DNA sequences: Computer simulations
- **Authors:** David Posada and Keith A. Crandall
- **Journal:** *Proceedings of the National Academy of Sciences USA* (2001), 98(24): 13757–13762
- **DOI:** [10.1073/pnas.241136198](https://doi.org/10.1073/pnas.241136198)
- **PMID:** 11717435

---

## Provenance Classification
- **Data Status:** **Faithfully Regenerated Synthetic Coalescent Suite**
- **Original Data Availability:** 
  The original 3,200 simulated alignments generated in 2001 were not deposited in an online public data repository (Dryad, Zenodo, and modern data archiving mandates did not exist at the time of publication).
- **Regeneration Protocol:** 
  The datasets were regenerated from scratch strictly adhering to the procedural and parametric specifications described in Posada & Crandall (2001):
  1. **Genealogical Simulation:** Coalescent simulations with crossing-over using Hudson's ancestral recombination graph model (`make_tree` / `ms`).
  2. **Sequence Evolution:** Finite-sites Markov substitution using `Seq-Gen` under the F84/HKY85 nucleotide substitution model with uniform base frequencies and transition/transversion ratio $\kappa = 2.0$.
  3. **Site-to-Site Rate Heterogeneity:** Modeled using a discrete 4-category Gamma distribution with shape parameter $\alpha \in \{\infty, 2.0, 0.5, 0.05\}$.

---

## Experimental Grid (3,200 Alignments Total)

All alignments comprise **$N = 10$ taxa** and **$L = 1{,}000$ nucleotides**, with **100 independent replicates** per scenario:

### 1. Recombination Power Scenarios (20 Scenarios, 2,000 Alignments)
- **Divergence parameter:** $\theta = 4N\mu \in \{10, 50, 100, 200\}$
- **Recombination intensity:** $\rho = 4Nr \in \{0, 1, 4, 16, 64\}$
- **Site rates:** Uniform ($\alpha = \infty$)
- **Coalescent partitions:** Ranges from 1.0 ($\rho = 0$) to a mean of $>160$ physical partitions per alignment ($\rho = 64$).

### 2. Rate Heterogeneity False Positive Scenarios (12 Scenarios, 1,200 Alignments)
- **Divergence parameter:** $\theta = 4N\mu \in \{10, 50, 100, 200\}$
- **Recombination intensity:** $\rho = 0$ (strictly non-recombinant null)
- **Gamma shape parameter:** $\alpha \in \{2.0\text{ (mild)}, 0.5\text{ (moderate)}, 0.05\text{ (severe)}\}$

---

## Methodological Comparisons
- **Direct Head-to-Head:** Both \RhizAeon (standard production engine with analytical crossover validation gate and profile likelihood polisher) and 3SEQ (exhaustive triplet hypergeometric test) are evaluated on the exact same 3,200 synthetic alignments.
- **Historical Benchmark Digits:** Historical method results (Maynard Smith Homoplasy Test, RETICULATE, MAXCHI, GENECONV, and RDP) were digitized directly from Posada & Crandall (2001, Figure 1).

---

## Associated Files
- [`generator.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/generator.py): Parametric coalescent and finite-sites sequence generator.
- [`scenarios.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/scenarios.py): Grid definitions matching Posada 2001.
- [`run_posada_benchmark.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/run_posada_benchmark.py): Parallel benchmarking harness for RhizAeon and 3SEQ.
- [`posada_raw_results.csv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/posada_raw_results.csv): Complete trial-by-trial results across all 3,200 runs.
- [`posada_summary.csv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/posada_summary.csv): Scenario-level summary statistics.
- [`plot_posada_benchmark.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/plot_posada_benchmark.py): Script generating Figure 1 replication plots.
- [`fig_posada_2001_replication.pdf`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/01_posada_2001/fig_posada_2001_replication.pdf): Multi-panel replication plot.
