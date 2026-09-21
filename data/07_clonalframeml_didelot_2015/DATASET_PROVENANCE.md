# Dataset Provenance: Didelot et al. (2015) ClonalFrameML & SimBac Bacterial Benchmark

## Canonical Citations
> **Didelot, X., & Wilson, D. J. (2015).**  
> *ClonalFrameML: Efficient Inference of Recombination in Whole Bacterial Genomes.*  
> **PLoS Computational Biology**, 11(2): e1004041.  
> DOI: [10.1371/journal.pcbi.1004041](https://doi.org/10.1371/journal.pcbi.1004041)  
> PMCID: [PMC4334968](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4334968/)

> **Didelot, X., & Falush, D. (2007).**  
> *Inference of bacterial microevolution using multilocus sequence data.*  
> **Genetics**, 175(3): 1251–1266.  
> DOI: [10.1534/genetics.106.063305](https://doi.org/10.1534/genetics.106.063305)  
> PMCID: [PMC1840082](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1840082/)

> **Brown, T., Didelot, X., Wilson, D. J., & Eyre, D. W. (2016).**  
> *SimBac: simulation of whole bacterial genomes with homologous recombination.*  
> **Microbial Genomics**, 2(1): e000044.  
> DOI: [10.1099/mgen.0.000044](https://doi.org/10.1099/mgen.0.000044)  
> PMCID: [PMC5320579](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5320579/)

---

## Dataset Architecture & Data Sources

This benchmark evaluates bacterial homologous gene conversion and microevolution across two complementary components:

### Component A: Original Empirical Whole Bacterial Genomes (Figshare)
- **Source:** Figshare public dataset DOI [10.6084/m9.figshare.19626912.v1](https://doi.org/10.6084/m9.figshare.19626912.v1), package `cfml.tgz` (97.8 MB).
- **Taxa & Chromosome:** 110 whole genomes of *Staphylococcus aureus* representing carriage and clinical reference isolates mapped to the MRSA252 reference genome ($L = 2{,}902{,}619$ nt; 106,480 variable single-nucleotide sites).
- **Known Biological Events:**
  - Sequence Type 34 (ST 34): 244-kb (inferred 231-kb) chromosomal replacement from an ST 10 donor spanning the origin of replication.
  - Sequence Type 239 (ST 239): 635-kb (inferred 555-kb) chromosomal replacement from an ST 30 donor into an ST 8 background.
  - Sequence Type 582 (ST 582): 310-kb novel chromosomal replacement spanning coordinates 845–1,155 kb (~1 Mb from origin).
- **Files Archived in `og_data/`:**
  - `Saureus.fasta`: 110 complete 2.9-Mb aligned sequences.
  - `Saureus.phyML.newick`: ML clonal genealogy estimated from core sites.
  - `Saureus.non-core-sites.txt`: 279,726 non-core coordinates excluded from phylogenetic core.
  - `core-sites.txt`: 2,622,893 core genome coordinates.
  - `cfml_results_2.R`: Authors' analysis and visualization script.

### Component B: Parametric Bacterial Microevolution Suite (ClonalFrame / SimBac)
- **Framework:** Coalescent simulation with point mutations ($\theta = 0.01$ per site) and Poisson homologous gene conversion imports at rate $R = (r/m \cdot \theta) / (\delta \cdot \nu)$ with mean tract length $\delta$ and sequence divergence $\nu$.
- **Experimental Design:** 22 parameter conditions (260 bacterial alignments, $L = 25{,}000$ nt):
  1. **Recombination Intensity ($r/m$):** $\{0.0, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0\}$ ($N=20, \delta=1000$ bp, $\nu=0.02$). Condition $r/m = 0.0$ serves as the clonal null negative control.
  2. **Tract Length ($\delta$):** $\{250, 500, 1000, 2500, 5000\}$ bp ($r/m=1.0, \nu=0.02, N=20$).
  3. **Import Divergence ($\nu$):** $\{0.005, 0.01, 0.03, 0.08\}$ ($r/m=1.0, \delta=1000$ bp, N=20).
  4. **Taxonomic Cohort Scaling ($N$):** $\{10, 20, 50, 100\}$ isolates ($r/m=1.0, \delta=1000$ bp, $\nu=0.02$).

---

## Standardized RhizAeon Execution Workflow
All alignments were evaluated under the standardized calibrated RhizAeon pipeline without parameter retuning:
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
