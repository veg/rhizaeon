# AGENT.MD: Autonomous Reproduction Protocol for the RhizAeon Benchmark Suite

**Authoritative, Deterministic Guide for AI Agents and Automated Pipelines to Reproduce all 10 Empirical and Simulated Recombination Benchmarks and Web Application from Scratch**

---

## 1. Executive Summary & Foundational Principles

This document provides complete, deterministic instructions for an autonomous AI agent (e.g., Antigravity, Claude, ChatGPT, Cursor) or automated CI/CD pipeline to execute, calibrate, audit, and reproduce the **10 curated empirical and simulated recombination benchmarks** (24,000+ genomes and alignments, spanning 1997–2026 CE) and build the complete publication-grade web application from scratch.

### Foundational Invariants for Autonomous Agents
1. **100% Tree-Free Continuous Coordinate Mechanics**:
   RhizAeon operates directly on prefix distance tensors $\mathcal{D}(s, i, j)$ embedded in continuous sequence manifolds. **Zero guide trees, zero bifurcating topologies, and zero tree-space branch-swapping moves are ever supplied, inferred, or traversed.**
2. **Authentic Author Alignment Preservation**:
   Whenever original author-deposited alignments exist (RDP5 empirical suites, RIVET pandemic trios, RecombinHunt 10,500 whole genomes, *S. aureus* 2.9-Mb pangenomes, H5N1 surveillance feeds), the pipeline strictly evaluates the original deposited data without gap alteration or simplex imputation.
3. **Strict Generative Replications**:
   Where historical benchmark alignments were not deposited at publication time (Posada 2001, PhiPack 2006, Olabode 2022), generative engines (`ms`, `Seq-Gen`, `msprime`, Wright-Fisher SPR) faithfully reproduce the published simulation protocols.
4. **Deterministic Reproducibility & Sub-Second Execution**:
   Every RhizAeon screen executes in 2.0 milliseconds to 10.37 seconds on commodity hardware, bypassing the combinatorial triplet growth ($\mathcal{O}(N^3)$) and tree partition searches ($\mathcal{O}(N^4)$) required by traditional tools.

---

## 2. System Requirements & Environment Setup

### 2.1 Hardware Architecture
- **Platforms**: Apple Silicon (M1/M2/M3/M4) or Linux/macOS x86_64 / ARM64.
- **RAM**: Minimum 8 GB available RAM (16 GB recommended for whole-chromosome cohorts like *S. aureus* [2.9 Mb]).
- **Storage**: ~2 GB free space for datasets, archives, diagnostics, figures, and web dossiers.

### 2.2 Python Dependencies
The execution environment requires Python $\ge 3.9$:
```bash
pip install numpy scipy pandas matplotlib biopython torch
```

### 2.3 Setting PYTHONPATH
RhizAeon must be in the Python module search path:
```bash
pip install -e .   # or: export PYTHONPATH="rhizaeon:$PYTHONPATH"
```

---

## 3. Directory Layout & Data Topology

```
rhizaeon_bench/
├── index.html                     # Master portal homepage with interactive search, filters, SVG scaling chart
├── benchmarks_master.json         # Master database: complete structured records for all 10 studies
├── AGENT.MD                       # Autonomous reproduction protocol & test harness
├── build_portal.py                # Self-contained portal builder & HTML generation script
├── study_curations.py             # Authoritative biological narratives & concordance taxonomy
├── studies/                       # 10 curated empirical & simulation benchmark dossiers
│   ├── 01_posada_2001/index.html
│   ├── ...
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

## 4. Benchmark Reproduction Procedures

To execute any individual benchmark study, unpack the corresponding data package and invoke the benchmark runner:

```bash
# Example: Replicating the canonical RDP5 empirical benchmark suite
cd data/04_rdp5_martin_2021
tar -xzf 04_rdp5_martin_2021_reproducibility.tar.gz
python3 run_rdp5_benchmark.py
python3 plot_rdp5_benchmark.py
```

### Rebuilding the Portal
After generating results or modifying study narratives:
```bash
python3 build_portal.py
```
This regenerates `index.html` and all 10 dossiers in `studies/<study_id>/index.html`.
