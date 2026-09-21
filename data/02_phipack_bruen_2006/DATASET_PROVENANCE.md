# Dataset Provenance: Bruen, Philippe, and Bryant (PhiPack 2006 Genetics)

## Reference Publication
- **Title:** A Simple and Robust Statistical Test for Detecting the Presence of Recombination
- **Authors:** Trevor C. Bruen, Hervé Philippe, and David Bryant
- **Journal:** *Genetics* (2006), 172(4): 2665–2681
- **DOI:** [10.1534/genetics.105.048975](https://doi.org/10.1534/genetics.105.048975)
- **PMID:** 16489234
- **PMCID:** PMC1456399

---

## Provenance Classification
- **Data Status:** **Faithfully Regenerated Synthetic Coalescent Suite**
- **Original Data Availability:**
  The synthetic datasets evaluated in Bruen et al. (2006) were generated using a modified version of *Treevolve* (Grassly et al. 1999) and Hudson's coalescent simulator. Alignments were not archived in public repositories in 2006.
- **Regeneration Protocol:**
  Alignments are regenerated strictly following the parametric and demographic specifications established in Bruen et al. (2006):
  1. **Genealogical Coalescent Modeling:** Neutral ancestral recombination graph (ARG) simulations under Hudson's formulation across constant-size populations and exponential population expansions ($N(t) = N_0 e^{-\beta t}$).
  2. **Mutational Parameterization:** Population mutation parameter $\theta = 4N\mu \in \{5, 10, 20\}$ per 1,000 bp locus.
  3. **Recombination Intensity:** $\rho = 4Nr \in \{0, 1, 4, 16, 64\}$.
  4. **Sample Sizes:** $n \in \{10, 50\}$ taxa; sequence length $L = 1{,}000$ bp.
  5. **Confounding Evolutionary Models ($\rho = 0$ Null):**
     - Exponential population growth producing star-like genealogical topologies.
     - Spatial substitution rate autocorrelation and Gamma rate heterogeneity ($\alpha \in \{0.5, 2.0\}$).

---

## Experimental Grid (2,300 Alignments Total)

### 1. Recombination Power Scenarios (1,200 Alignments)
- Sample sizes: $n \in \{10, 50\}$ taxa
- Mutation levels: $\theta \in \{5, 10, 20\}$
- Recombination intensities: $\rho \in \{1, 4, 16, 64\}$
- Sequence length: $L = 1{,}000$ nt
- Replicates: 50 independent replicates per configuration ($2 \times 3 \times 4 \times 50 = 1{,}200$ alignments)

### 2. Rate Variation & Demographic False Positive Controls (1,100 Alignments, $\rho = 0$)
- Constant-size neutral null: $n \in \{10, 50\}$, $\theta \in \{5, 10, 20\}$ (300 alignments)
- Exponential growth null: $n \in \{10, 50\}$, $\theta \in \{10, 20\}$, growth rate $\beta \in \{5.0, 20.0\}$ (400 alignments)
- Rate heterogeneity null: $n \in \{10, 50\}$, $\theta \in \{10, 20\}$, Gamma shape $\alpha \in \{0.5, 2.0\}$ (400 alignments)

---

## Evaluated Methods
- **RhizAeon:** Standard default production workflow (prefix distance tensor + scale-aware bilateral RP-FDA screen + analytical crossover validation gate $p < 0.005$ + profile likelihood polisher).
- **PHI ($\Phi_w$):** Pairwise Homoplasy Index test (Bruen et al. 2006, executed via the official `PhiPack` binary).
- **Max $\chi^2$:** Maximum $\chi^2$ test (Maynard Smith 1992, executed via `PhiPack -o`).
- **NSS:** Neighbor Similarity Score (Jakobsen & Easteal 1996, executed via `PhiPack -o`).
- **3SEQ:** Exhaustive non-parametric triplet mosaic test (Boni et al. 2007).

---

## Associated Files
- [`generator.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_phipack_bruen_2006/generator.py): Parametric coalescent simulator with exponential growth and rate heterogeneity.
- [`run_phipack_benchmark.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_phipack_bruen_2006/run_phipack_benchmark.py): Parallel benchmarking harness executing RhizAeon, PHI, MaxChi, NSS, and 3SEQ.
- [`phipack_raw_results.csv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_phipack_bruen_2006/phipack_raw_results.csv): Complete trial-by-trial results across all 2,300 runs.
- [`phipack_summary.csv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_phipack_bruen_2006/phipack_summary.csv): Scenario-level summary metrics.
- [`plot_phipack_benchmark.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_phipack_bruen_2006/plot_phipack_benchmark.py): Production visualization script.
- [`fig_phipack_replication.pdf`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_phipack_bruen_2006/fig_phipack_replication.pdf): Publication-ready multi-panel figure.
