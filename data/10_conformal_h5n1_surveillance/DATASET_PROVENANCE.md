# Dataset Provenance: Conformal Predictive Outbreak Surveillance (7,054 Complete H5N1 Genomes)

## Canonical Citations & Conceptual Framework
> **Vovk, V., Gammerman, A., & Shafer, G. (2005).**  
> *Algorithmic Learning in a Random World.*  
> Springer Science & Business Media.  
> DOI: [10.1007/b106715](https://doi.org/10.1007/b106715)

> **Kosakovsky Pond, S. L., et al. (2026).**  
> *Clean Breaks: Tree-Free Discovery of Mosaic Genomes on Continuous Sequence Manifolds.*  
> In preparation.

---

## Dataset Architecture & Empirical Cohort Provenance

This study evaluates distribution-free conformal prediction and streaming reassortment surveillance across the complete global H5N1 clade 2.3.4.4b panzootic outbreak:

1. **Empirical Panzootic Archive:**
   - 7,054 complete eight-segment Influenza A virus genomes (56,432 segments, 89.9 million nucleotides) collected globally between 1996 and 2024 from NCBI GenBank and GISAID.
   - Spans twelve host categories: dairy cattle ($N = 866$), marine mammals ($N = 43$), seabirds & shorebirds ($N = 303$), domestic poultry, domestic cats ($N = 60$), red foxes ($N = 96$), and human farmworkers ($N = 2$).

2. **Conformal Calibration Design:**
   - Calibration set $N_{\mathrm{calib}} = 400$ diversity genomes selected via stratified furthest-point medoid sampling across all twelve host interfaces.
   - Non-conformity score: Normalized pairwise Hamming divergence on metric manifolds against nearest established lineage medoids:
     $$s_q^{(s)} = \min_{m \in \mathcal{M}_s} \mathcal{D}_s(q, m)$$
   - Segment p-values evaluated on the exact finite-sample conformal lattice:
     $$p^{(s)}(q) = \frac{1 + \sum_{j=1}^{N_{\mathrm{calib}}} \mathbf{1}\left(s_j^{(s)} \ge s_q^{(s)}\right)}{N_{\mathrm{calib}} + 1}$$
   - Whole-genome intersection hypothesis tested via the Simes omnibus test:
     $$p_{\mathrm{simes}}(q) = \min_{k=1,\dots,8} \left\{ \frac{8}{k} p_{(k)}(q) \right\}$$

3. **Key Findings:**
   - Exact mathematical error calibration: At nominal $\alpha = 0.05$, exactly 261 genomes (3.70\%) reject conformity (within the $\le 5\%$ error budget). At $\alpha = 0.01$, exactly 49 genomes (0.69\%) reject.
   - Host novelty hotspots: Marine mammals exhibit 16.28\% novelty, seabirds 7.26\%, domestic cats 6.67\%, while dairy cattle show 95.61\% conformity and human farmworker cases show intact acquisition of the bovine B3.13 cassette ($p > 0.48$).
   - Dynamic Medoid Recruitment suppresses alarm fatigue by 44.4\% (from 261 down to 145 flagged genomes), recruiting 206 exemplar medoids across segments.
   - Wall-clock latency: 0.095 seconds across all 56,432 segments (0.013 ms per complete genome, >70,000 genomes/sec).

---

## Standardized Execution Pipeline
Evaluated using `rhizaeon.conformal.ConformalPredictor` and `PrefixDistanceEngine` across the 8-segment genome.
