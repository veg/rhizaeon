# Dataset Provenance: Alfonsi et al. (2024) RecombinHunt Benchmark

## Canonical Citation
> **Alfonsi, E., Cereda, M., & Pelizzola, M. (2024).**  
> *RecombinHunt: a compressed-sensing-based method for identifying recombinant viral genomes and breakpoint locations.*  
> **Nature Communications**, 15(1), 3717.  
> DOI: [10.1038/s41467-024-47464-5](https://doi.org/10.1038/s41467-024-47464-5)  
> PMCID: [PMC11065798](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11065798/)

---

## Dataset Architecture & Original Benchmark Alignments

This benchmark evaluates 10,500 synthetic SARS-CoV-2 complete genomes ($L = 29{,}903$ nt) distributed directly in official Supplementary Data 1 (`supplementary_data_1.txt`, 12.7 MB) alongside empirical Mpox genomes in Supplementary Data 2 (`supplementary_data_2_mpox.txt`):

1. **0-Breakpoint Clonal Nulls (3,500 genomes):**
   - Pure vertical mutation models evaluating Type I error control under increasing mutational noise.
   - Seven noise tiers: 1, 3, 5, 10, 15, 20, and 30 random private mutations per genome (500 replicates each).
   - Backbone: Omicron BA.2.

2. **1-Breakpoint Single Mosaics (3,500 genomes):**
   - Single-crossover recombinant genomes between Omicron BA.2 and Delta AY.45.
   - Seven noise tiers: 0, 3, 5, 10, 15, 20, and 30 random private mutations per genome (500 replicates each).

3. **2-Breakpoint Cassette Mosaics (3,500 genomes):**
   - Double-crossover gene conversion cassettes (Omicron BA.2 backbone with an inserted Delta AY.45 genomic tract).
   - Seven noise tiers: 0, 3, 5, 10, 15, 20, and 30 random private mutations per genome (500 replicates each).

4. **Empirical Mpox Cohort (Supplementary Data 2):**
   - Authentic Orthopoxvirus monkeypox genomes.

---

## Benchmarked Methods & Baseline Sources
1. **RhizAeon:** Evaluated locally on all 10,500 genomes using the calibrated standard workflow:
   `PrefixDistanceEngine(mat, codon_aligned=False, compute_transitions=True)`
   `run_recursive_partition_fda_screen(engine, taxa_names=taxa, min_z=1.8, min_pir=0.06, crossover_validation=True, crossover_p_threshold=0.005, min_informative_sites=3, polish_ml=True)`
2. **RecombinHunt:** Digitized directly from Alfonsi et al. (2024) Supplementary Data 1 and Tables 1a, 1b.
