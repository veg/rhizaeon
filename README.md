# RhizAeon Benchmark Compendium & Reproducibility Portal

This repository hosts the static, publication-grade web application and complete reproducibility compendium documenting the **10 curated empirical and simulated benchmark cohorts** (24,000+ genomes and alignments, 1997–2026) evaluated in the **RhizAeon** manuscript:

> **"Clean Breaks: Real-Time Recombination Detection, Exact Coordinate Geometry, and the Demarcation of Viral Chimeras"**  
> *Sergei L. Kosakovsky Pond, ... , Darren P. Martin.*

---

## 1. Live Interactive Web Compendium

Explore the full benchmark results, multi-panel diagnostic figures, interactive charts, and downloadable data archives online:

* **Master Portal:** [https://veg.github.io/rhizaeon/](https://veg.github.io/rhizaeon/)
* **RecombinHunt 10,500 Full-Genome Dossier:** [`studies/06_recombinhunt_alfonsi_2024/index.html`](https://veg.github.io/rhizaeon/studies/06_recombinhunt_alfonsi_2024/index.html)
* **Panzootic H5N1 7,054-Genome Conformal Surveillance Dossier:** [`studies/10_conformal_h5n1_surveillance/index.html`](https://veg.github.io/rhizaeon/studies/10_conformal_h5n1_surveillance/index.html)
* **Canonical RDP5 Empirical Multi-Virus Dossier:** [`studies/04_rdp5_martin_2021/index.html`](https://veg.github.io/rhizaeon/studies/04_rdp5_martin_2021/index.html)

---

## 2. Compendium Highlights & Key Metrics

* **10 Curated Empirical & Simulated Cohorts**: Incorporates authentic author-deposited alignments (RDP5 empirical suites, RIVET pandemic trios, RecombinHunt 10,500 whole genomes, *S. aureus* 2.9-Mb pangenomes, H5N1 surveillance feeds) alongside strict generative coalescent replications (Posada 2001, PhiPack 2006, Olabode 2022).
* **100% Tree-Free Continuous Coordinate Manifolds**: Completely bypasses guide tree reconstruction, branch swapping, and combinatorial triplet enumeration, replacing discrete bifurcations with continuous sequence geometry and prefix distance tensors.
* **Sub-Second to Sub-Minute Execution**: Computes in 2.0 milliseconds to 10.37 seconds on commodity hardware across viral and bacterial cohorts, yielding **5.6&times; to 142.6&times; speedups over RDP5**, **up to 798&times; over RDP4**, and **>15,000&times; over phylogenetic partitioning (GARD)**.
* **Zero-Nucleotide Exact Physical Coordinate Recovery**: The Data-Driven Maximum Likelihood Polisher resolves the uninformative plateau $[x_{\mathrm{left}}, x_{\mathrm{right}}]$ between flanking SNPs, achieving exact 0-nt concordance with Sanger-sequenced breakpoints in HIV-1 KAL153, Potato Virus Y, and Begomovirus TYLCV.
* **Strict Empirical Error Control**: Analytical crossover validation gates strictly bound false positive rates to $\le 2.0\%$ under extreme substitution rate heterotachy ($\alpha = 0.05$) and explosive star-like demographic expansions ($\beta = 20$), regimes where historical compatibility statistics and homoplasy tests collapsed with 28% to 90% false alarms.
* **100% Verified Literature Provenance**: Every study maintains canonical DOIs, CrossRef-verified citations, author-deposited data manifests, and standalone replication scripts.

---

## 3. Planetary-Scale Surveillance Grand Challenges

Beyond small-to-medium alignments, RhizAeon tackles the real-world operational challenges of streaming genomic surveillance:

### Challenge A: The RecombinHunt 10,500 Full-Genome Grand Challenge
* **Dataset:** 10,500 full-length 29,903-nucleotide SARS-CoV-2 whole genomes from Alfonsi et al. (2024) *Nat. Commun.* across 21 parameter cells (0 to 30 private noise mutations).
* **Head-to-Head Comparison:** Evaluated against RecombinHunt and GARD.
* **Results:**
  - Evaluates the entire 10,500-genome compendium in **21.2 seconds total** (**2.02 ms per genome; 495 genomes/second throughput**).
  - **0.00% False Positive Rate** across 3,500 non-recombinant control genomes (where RecombinHunt false alarms escalated to 8.8%).
  - **100.0% Detection Sensitivity** and **100.0% exact flanking interval placement** on single-crossover mosaics (1-BP).
  - **99.8% Sensitivity** and **94.2% exact placement** on double-crossover cassettes (2-BP) under 30 noise substitutions (vs 78.8% for RecombinHunt).
* **Dossier:** [`studies/06_recombinhunt_alfonsi_2024/index.html`](https://veg.github.io/rhizaeon/studies/06_recombinhunt_alfonsi_2024/index.html)

### Challenge B: The Panzootic Avian Flu (H5N1) 7,054 Complete Genome Conformal Screen
* **Dataset:** 7,054 complete 8-segment H5N1 genomes (89.9 megabases, 56,432 individual segments) harvested from NCBI and GISAID surveillance feeds.
* **Operational Barrier:** Replacing multi-segment tree building with real-time anomaly triage.
* **Results:**
  - **Real-Time Streaming Triage**: Evaluates complete 8-segment constellations in **0.94 seconds per genome**.
  - **Conformal Guarantees**: Achieves **99.4% empirical coverage** at nominal $\alpha = 0.01$ (bounded 0.6% error rate) under distribution-free non-conformity scoring on metric Grassmannian spaces.
  - **Epizootic Constellation Tracking**: Automatically isolates the B3.13 dairy cattle epizootic genotype in 15/15 bovine samples (100.0% prevalence) and identifies novel inter-clade reassortants in domestic cats, red foxes, and marine mammals.
  - **Dynamic Medoid Recruitment**: Continuously updates reference medoids to absorb expanding transmission waves and prevent surveillance alarm fatigue.
* **Dossier:** [`studies/10_conformal_h5n1_surveillance/index.html`](https://veg.github.io/rhizaeon/studies/10_conformal_h5n1_surveillance/index.html)

---

## 4. Master Benchmark Results Table

| **#** | **Cohort & Primary Reference** | **Category** | **Taxa ($N$)** | **Length ($L$)** | **Scale** | **Baseline Runtime** | **RhizAeon Latency** | **Speedup** | **Error Control / Power** | **Dossier** |
| :---: | :--- | :--- | :---: | :---: | :---: | :--- | :--- | :---: | :--- | :---: |
| **01** | **Posada & Crandall (2001)** *PNAS* | Coalescent | 10 | 1,000 nt | 3,200 alns | 39.52 ms (3SEQ) | **8.92 ms** | **4.4&times;** | 2.0% FPR under $\alpha=0.05$; 100% power | [View &rarr;](studies/01_posada_2001/index.html) |
| **02** | **Bruen et al. (PhiPack 2006)** *Genetics* | Coalescent | 10–50 | 1,000 nt | 2,300 alns | 340.7 ms (PhiPack) | **38.8 ms** | **8.8&times;** | 1.0% FPR under $\beta=20$; 96% power | [View &rarr;](studies/02_phipack_bruen_2006/index.html) |
| **03** | **Olabode et al. (2022)** *PNAS* | Phylodynamics | 16–200 | 9,000 nt | 360 alns | 22 min (RDP5) / 138 min (RDP4) | **10.37 s** | **127&times; / 798&times;** | 100% sensitivity; 28.5 nt median MAE | [View &rarr;](studies/03_olabode_2022/index.html) |
| **04** | **Martin et al. (RDP5 2021)** *MBE* | Multi-Virus | 6–274 | 2,989–9,953 nt | 6 cohorts | 1.9 s–1.02 h (RDP5) | **0.208 s–10.85 min** | **5.6&times;–142.6&times;** | 0-nt exact Sanger breakpoint concordance | [View &rarr;](studies/04_rdp5_martin_2021/index.html) |
| **05** | **RIVET (Smith et al. 2023)** *Bioinformatics* | SARS-CoV-2 | 3 (trios) | 29,903 nt | 2,381 genomes | 49.3 ms (3SEQ) | **5.4 ms** | **9.1&times;** | 24.5% power at $d=10$ (3SEQ: 0.0%); 0.12% FPR | [View &rarr;](studies/05_rivet_sars_cov_2/index.html) |
| **06** | **RecombinHunt (Alfonsi 2024)** *Nat. Commun.* | SARS-CoV-2 | 10,500 | 29,903 nt | 10,500 genomes | ~8 min (RecombinHunt) | **21.2 s (total)** | **22.5&times;** | 0.00% FPR (3,500 nulls); 100% 1-BP accuracy | [View &rarr;](studies/06_recombinhunt_alfonsi_2024/index.html) |
| **07** | **ClonalFrameML (Didelot 2015)** *PLoS Comp. Biol.* | Bacteria | 110 + 20 | 2.9 Mb / 25 kb | 110 + 260 sims | 15–60 min (CFML) | **8.4 s (110 S. aureus)** | **>150&times;** | 93.4% within $\pm 50$ nt for micro-imports | [View &rarr;](studies/07_clonalframeml_didelot_2015/index.html) |
| **08** | **Ghost Introgression Suite** *Current Study* | Deep Introgression | 5–100 | 3,000 nt | 370 alns | Hours ($\mathcal{O}(N^3)$ search) | **608 ms (N=100)** | **266&times;** | 95–100% power under 5% donor drift | [View &rarr;](studies/08_ghost_introgression_suite/index.html) |
| **09** | **Grand 1,000 Codon Suite** *Current Study* | Codon Selection | 16 | 1,800 nt | 1,000 alns | 109.1 ms (3SEQ) | **11.5 ms (RP-FDA)** | **9.5&times;** | 0.0% FPR under intense positive selection ($\omega=5$) | [View &rarr;](studies/09_selection_vs_recombination/index.html) |
| **10** | **Conformal H5N1 Surveillance** *Current Study* | Surveillance | 7,054 | 8 segments | 7,054 genomes | Hours (ML trees) | **0.94 s / genome** | **Real-time** | 99.4% coverage; 0.6% conformal error rate | [View &rarr;](studies/10_conformal_h5n1_surveillance/index.html) |

---

## 5. Directory Structure & File Map

```
rhizaeon_bench/
├── index.html                     # Master portal homepage with interactive search, filters, SVG scaling chart
├── benchmarks_master.json         # Master database: complete structured records for all 10 studies
├── AGENT.MD                       # Comprehensive autonomous agent reproduction protocol & test harness
├── build_portal.py                # Self-contained portal builder & HTML generation script
├── study_curations.py             # Authoritative biological narratives & concordance taxonomy
├── studies/                       # 10 curated empirical & simulation benchmark dossiers
│   ├── 01_posada_2001/index.html
│   ├── 02_phipack_bruen_2006/index.html
│   ├── 03_olabode_2022/index.html
│   ├── 04_rdp5_martin_2021/index.html
│   ├── 05_rivet_sars_cov_2/index.html
│   ├── 06_recombinhunt_alfonsi_2024/index.html
│   ├── 07_clonalframeml_didelot_2015/index.html
│   ├── 08_ghost_introgression_suite/index.html
│   ├── 09_selection_vs_recombination/index.html
│   └── 10_conformal_h5n1_surveillance/index.html
├── data/                          # Complete primary data & reproducibility artifacts (.tar.gz)
│   ├── 01_posada_2001/01_posada_2001_reproducibility.tar.gz
│   ├── ...
│   └── 10_conformal_h5n1_surveillance/10_conformal_h5n1_surveillance_reproducibility.tar.gz
├── assets/
│   ├── css/style.css              # BRC-Analytics design system & responsive layouts
│   ├── js/main.js                 # Client-side filtering, interactive SVG scaling chart, lightbox
│   ├── figures/                   # Publication-grade vector PDF and high-res PNG figures
│   └── img/                       # Walkthrough figures and diagrams
```

---

## 6. Installation and CLI Usage

Install RhizAeon via `pip`:

```bash
pip install rhizaeon
```

### Basic Command-Line Invocations

```bash
# 1. Fast screening of an aligned viral cohort
rhizaeon screen --input alignment.fasta --output results.csv

# 2. Maximum likelihood profile polishing of candidate breakpoints
rhizaeon polish --input alignment.fasta --breakpoints results.csv --output polished.csv

# 3. Conformal multi-segment reassortment triage
rhizaeon conformal --segments segment_manifest.tsv --alpha 0.01 --output triage_report.csv
```

---

## 7. Citation

If you use RhizAeon or this benchmark compendium in your research, please cite:

```bibtex
@article{rhizaeon2026,
  title={Clean Breaks: Real-Time Recombination Detection, Exact Coordinate Geometry, and the Demarcation of Viral Chimeras},
  author={Kosakovsky Pond, Sergei L. and Martin, Darren P. and others},
  journal={bioRxiv},
  year={2026},
  publisher={Cold Spring Harbor Laboratory}
}
```
