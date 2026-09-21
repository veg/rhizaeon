# Dataset Provenance: Selection vs Recombination (Codon Spectrum Benchmark)

## Canonical Citations & Conceptual Framework
> **Anisimova, M., Nielsen, R., & Yang, Z. (2003).**  
> *Effect of recombination on the accuracy of the likelihood method for detecting positive selection at amino acid sites.*  
> **Genetics**, 164(3), 1229–1236.  
> DOI: [10.1093/genetics/164.3.1229](https://doi.org/10.1093/genetics/164.3.1229)  
> PMCID: [PMC1462615](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1462615/)

> **Kosakovsky Pond, S. L., Posada, D., Gravenor, M. B., Woelk, C. H., & Frost, S. D. (2006).**  
> *GARD: a genetic algorithm for recombination detection.*  
> **Bioinformatics**, 22(24), 3096–3098.  
> DOI: [10.1093/bioinformatics/btl474](https://doi.org/10.1093/bioinformatics/btl474)

> **Spielman, S. J., & Wilke, C. O. (2015).**  
> *Pyvolve: A flexible Python module for simulating sequences along phylogenies.*  
> **PLoS ONE**, 10(9), e0139047.  
> DOI: [10.1371/journal.pone.0139047](https://doi.org/10.1371/journal.pone.0139047)  
> PMCID: [PMC4575133](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4575133/)

---

## Benchmark Design & Experimental Conditions

A longstanding challenge in evolutionary bioinformatics is the confounding between horizontal genetic exchange and positive diversifying selection ($dN/dS = \omega > 1$). Episodic bursts of nonsynonymous substitutions create localized clusters of homoplasy that mimic crossover breakpoints in sliding-window and phylogenetic tests. Conversely, unmodeled recombination inflates false positive rates in codon models of positive selection (Anisimova et al. 2003).

This benchmark evaluates 1,000 mechanistic codon alignments ($N = 16$ taxa, $L = 600$ codons / 1,800 nt, 50 stochastic replicates per cell) across 20 evolutionary scenarios generated via Pyvolve:

1. **Selection & Rate Heterogeneity Null Controls (Dimension A, 200 Alignments):**
   - Pure positive selection ($\omega = 0.85$ and episodic bursts with $\omega > 1.0$) under strictly clonal vertical descent.
   - Spatial rate variation (Gamma $\alpha = 0.3$ and $\alpha = 0.8$).
   - Homogeneous neutral evolution baseline.
   - Evaluates Type I error suppression (avoiding false mosaic calls induced by selection clusters).

2. **Divergence Spectrum & Saturation (Dimension B, 250 Alignments):**
   - Tree scales spanning two orders of magnitude: $T \in [0.05, 4.00]$ substitutions per site.

3. **Recombination Event Structural Complexity (Dimension C, 200 Alignments):**
   - Single crossovers, double-crossover cassettes, triple mosaics, and quadruple mosaics.

4. **Tract Length Sensitivity (Dimension D, 200 Alignments):**
   - Tract lengths from micro-exchanges ($\ell = 30$ codons / 90 nt) to standard cassettes ($\ell = 300$ codons / 900 nt).

5. **Parental Donor Architecture (Dimension E, 150 Alignments):**
   - Close vs distant ghost donors and intra-clade vs inter-clade parentage.

---

## Standardized RhizAeon Execution Workflow
Evaluations executed under the calibrated standard RhizAeon pipeline:
```python
engine = PrefixDistanceEngine(mat, codon_aligned=True, compute_transitions=True)
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
