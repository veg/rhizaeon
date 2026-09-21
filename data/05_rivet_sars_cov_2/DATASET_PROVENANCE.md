# Dataset Provenance: RIVET SARS-CoV-2 Recombination Benchmark

## Reference Publications
- **Primary Tool Publication:**
  - **Title:** Tracking and curating putative SARS-CoV-2 recombinants with RIVET
  - **Authors:** Kyle Smith, Cheng Ye, and Yatish Turakhia
  - **Journal:** *Bioinformatics* (2023), 39(9): btad538
  - **DOI:** [10.1093/bioinformatics/btad538](https://doi.org/10.1093/bioinformatics/btad538)
  - **PMID:** 37651464
  - **PMCID:** PMC10493179
- **Algorithmic Engine & Benchmark Protocol Publication:**
  - **Title:** Pandemic-scale phylogenomics reveals the SARS-CoV-2 recombination landscape
  - **Authors:** Yatish Turakhia, Bryan Thornlow, Angie S. Hinrichs, Nicola De Maio, Landen Gozashti, Robert Lanfear, David Haussler, and Russell Corbett-Detig
  - **Journal:** *Nature* (2022), 609(7929): 994–997
  - **DOI:** [10.1038/s41586-022-05189-9](https://doi.org/10.1038/s41586-022-05189-9)
  - **PMID:** 35921869
  - **PMCID:** PMC9519458

---

## Provenance Classification
This benchmark incorporates a dual-component design combining original empirical pandemic data and faithfully regenerated parametric simulations:

### Component A: Original Empirical Pandemic Recombinants (OG Data)
- **Data Status:** **Original Empirical Dataset from Official Public Repository**
- **Data Source:** RIVET Public Database (`https://rivet.ucsd.edu/download_public_table`), downloaded directly from the active surveillance portal.
- **Dataset Composition:**
  - Total screened recombinant candidate trios: 6,432 trios spanning millions of public global SARS-CoV-2 genomes.
  - High-confidence quality-control filtered trios: **481 trios passing all rigorous automated quality checks (`QC Flags == PASS`)**.
  - Includes 62 designated Pango recombinant lineages (e.g., XBB, XFG, XFV, XFQ, etc.).
  - Each entry contains complete coordinates of informative segregating sites, donor/acceptor/recombinant lineage alleles across the 29,903 nt Wuhan-Hu-1 reference genome, parsimony score improvements, RIVET-inferred breakpoint intervals, and author-computed 3SEQ statistics.

### Component B: Faithfully Regenerated SARS-CoV-2 Parametric Simulations
- **Data Status:** **Faithfully Regenerated Synthetic Benchmark Suite**
- **Regeneration Protocol:**
  - Adheres strictly to the simulation protocol developed by Bryan Thornlow and Yatish Turakhia (`bpt26/recombination/makeRandomRecombinants.py`), used to benchmark RIPPLES and RIVET in *Nature* (2022) and *Bioinformatics* (2023).
  - Background sequence: Canonical SARS-CoV-2 Wuhan-Hu-1 reference genome ($L = 29{,}903$ nt).
  - Recombination events: 1-breakpoint and 2-breakpoint mosaic genomes with breakpoints chosen uniformly at random across the 29,903 nt genome.
  - Parental mutational divergence: $d \in \{10, 20, 30, 50\}$ nucleotide substitutions separating donor and acceptor clades.
  - Post-recombination mutations: $m \in \{0, 1, 2, 3\}$ private substitutions accumulated after recombination.
  - Non-recombinant negative control clades: Identical mutational depths ($d \in \{10, 20, 30, 50\}$, $m \in \{0, 1, 2, 3\}$) under pure vertical mutation without recombination, used to assess empirical false positive rates.
  - Replicates: $N = 50$ independent replicates per grid condition (800 recombinant alignments and 400 negative control alignments, totaling 1,200 simulated whole genomes).

---

## Methodological Comparisons
- **Evaluated Tools:**
  1. **RhizAeon:** Standard default production workflow (prefix distance tensor + scale-aware bilateral RP-FDA screen + analytical crossover validation gate $p < 0.005$ + profile likelihood polisher).
  2. **3SEQ:** Exhaustive non-parametric triplet hypergeometric mosaic test (standard reference method incorporated directly inside the RIVET platform).
  3. **RIVET / RIPPLES:** Mutation-Annotated Tree (MAT) partial parsimony placement engine.
- **Key Evaluated Metrics:**
  - Power / sensitivity as a function of parental divergence $d$ and post-recombination mutations $m$.
  - False discovery rate under non-recombinant vertical descent.
  - Breakpoint interval concordance on empirical pandemic trios and spatial error $|\hat{b} - b_{\text{true}}|$ on simulated genomes.
  - Computational throughput and per-genome runtime scaling across 30 kb viral chromosomes.

---

## Associated Files
- [`data/rivet_public_recombinants.tsv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/data/rivet_public_recombinants.tsv): Original empirical RIVET public database table.
- [`data/wuhan.ref.fa`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/data/wuhan.ref.fa): Canonical SARS-CoV-2 reference sequence.
- [`generator.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/generator.py): Faithful simulation generator replicating Thornlow & Turakhia's protocol.
- [`run_rivet_benchmark.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/run_rivet_benchmark.py): Unified benchmarking engine evaluating RhizAeon and 3SEQ on empirical trios and synthetic grids.
- [`rivet_raw_results.csv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/rivet_raw_results.csv): Comprehensive trial-by-trial results.
- [`rivet_summary.csv`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/rivet_summary.csv): Stratified performance summary table.
- [`plot_rivet_benchmark.py`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/plot_rivet_benchmark.py): Production visualization script.
- [`fig_rivet_benchmark.pdf`](file:///Users/sergei/Projects/TOGA_MEME/recombination/simulations/02_rivet_sars_cov_2/fig_rivet_benchmark.pdf): Publication-ready multi-panel figure.
