#!/usr/bin/env python3
"""
build_portal.py
Authoritative portal builder script that compiles the 10 canonical RhizAeon benchmarks
into the publication-grade web application in rhizaeon_bench.
"""

import os
import sys
import json
import shutil
import math
import html

# Ensure local imports work
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from study_curations import STUDY_CURATIONS

PORTAL_DIR = os.path.dirname(os.path.abspath(__file__))
STUDIES_DIR = os.path.join(PORTAL_DIR, "studies")
ASSETS_DIR = os.path.join(PORTAL_DIR, "assets")
DATA_DIR = os.path.join(PORTAL_DIR, "data")
FIGURES_DIR = os.path.join(ASSETS_DIR, "figures")

# Define scaling points for the interactive SVG plot (N * L complexity vs runtime ms)
SCALING_DATA = [
    {
        "id": "01_posada_2001",
        "name": "Posada & Crandall (2001) PNAS",
        "category": "Classical Coalescent",
        "complexity_label": "10 taxa × 1,000 nt = 10,000 nt",
        "log_complexity": 4.00,
        "rhiz_time": "8.92 ms",
        "log_rhiz_ms": 0.95,
        "hist_time": "39.52 ms",
        "log_hist_ms": 1.60,
        "hist_method": "3SEQ",
        "speedup": "4.4x vs 3SEQ"
    },
    {
        "id": "02_phipack_bruen_2006",
        "name": "Bruen et al. (2006) PhiPack",
        "category": "Classical Coalescent",
        "complexity_label": "50 taxa × 1,000 nt = 50,000 nt",
        "log_complexity": 4.70,
        "rhiz_time": "38.8 ms",
        "log_rhiz_ms": 1.59,
        "hist_time": "340.7 ms",
        "log_hist_ms": 2.53,
        "hist_method": "PhiPack",
        "speedup": "8.8x vs PhiPack"
    },
    {
        "id": "03_olabode_2022",
        "name": "Olabode et al. (2022) HIV-1",
        "category": "Continuous-Time & Phylodynamics",
        "complexity_label": "200 taxa × 9,000 nt = 1,800,000 nt",
        "log_complexity": 6.26,
        "rhiz_time": "10.37 s",
        "log_rhiz_ms": 4.02,
        "hist_time": "1,320 s (22 min)",
        "log_hist_ms": 6.12,
        "hist_method": "RDP5 (138 min RDP4)",
        "speedup": "127x vs RDP5 (798x vs RDP4)"
    },
    {
        "id": "04_rdp5_martin_2021",
        "name": "Martin et al. (2021) FMDV",
        "category": "Empirical Multi-Virus",
        "complexity_label": "91 taxa × 8,299 nt = 755,209 nt",
        "log_complexity": 5.88,
        "rhiz_time": "1.62 s",
        "log_rhiz_ms": 3.21,
        "hist_time": "231 s (3.85 min)",
        "log_hist_ms": 5.36,
        "hist_method": "RDP5 (15.1 min RDP4)",
        "speedup": "142.6x vs RDP5"
    },
    {
        "id": "05_rivet_sars_cov_2",
        "name": "RIVET SARS-CoV-2 (Smith 2023)",
        "category": "SARS-CoV-2 Pandemic",
        "complexity_label": "3 taxa × 29,903 nt = 89,709 nt",
        "log_complexity": 4.95,
        "rhiz_time": "5.4 ms",
        "log_rhiz_ms": 0.73,
        "hist_time": "49.3 ms",
        "log_hist_ms": 1.69,
        "hist_method": "3SEQ",
        "speedup": "9.1x vs 3SEQ"
    },
    {
        "id": "06_recombinhunt_alfonsi_2024",
        "name": "Alfonsi et al. (2024) RecombinHunt",
        "category": "SARS-CoV-2 Pandemic",
        "complexity_label": "10,500 taxa × 29,903 nt = 313,981,500 nt",
        "log_complexity": 8.50,
        "rhiz_time": "21.2 s (total batch)",
        "log_rhiz_ms": 4.33,
        "hist_time": "472.5 s (~8 min)",
        "log_hist_ms": 5.67,
        "hist_method": "RecombinHunt (>15,000x vs GARD)",
        "speedup": "22.5x vs RecombinHunt"
    },
    {
        "id": "07_clonalframeml_didelot_2015",
        "name": "Didelot et al. (2015) S. aureus",
        "category": "Bacterial Microevolution",
        "complexity_label": "110 taxa × 2,902,619 nt = 319,288,090 nt",
        "log_complexity": 8.50,
        "rhiz_time": "8.4 s (106k SNPs)",
        "log_rhiz_ms": 3.92,
        "hist_time": "1,800 s (30 min)",
        "log_hist_ms": 6.25,
        "hist_method": "ClonalFrameML",
        "speedup": ">150x vs ClonalFrameML"
    },
    {
        "id": "08_ghost_introgression_suite",
        "name": "Ghost Introgression Suite",
        "category": "Deep Phylogenomics & Introgression",
        "complexity_label": "100 taxa × 3,000 nt = 300,000 nt",
        "log_complexity": 5.48,
        "rhiz_time": "608.0 ms",
        "log_rhiz_ms": 2.78,
        "hist_time": "161,700 ms (161.7 s)",
        "log_hist_ms": 5.21,
        "hist_method": "O(N^3) Triplet Enumeration",
        "speedup": "266x vs Triplet Search"
    },
    {
        "id": "09_selection_vs_recombination",
        "name": "Selection vs Recombination Grand 1,000",
        "category": "Codon Selection & Architecture",
        "complexity_label": "16 taxa × 1,800 nt = 28,800 nt",
        "log_complexity": 4.46,
        "rhiz_time": "11.5 ms (RP-FDA)",
        "log_rhiz_ms": 1.06,
        "hist_time": "109.1 ms",
        "log_hist_ms": 2.04,
        "hist_method": "3SEQ",
        "speedup": "9.5x vs 3SEQ"
    },
    {
        "id": "10_conformal_h5n1_surveillance",
        "name": "Conformal H5N1 Surveillance",
        "category": "Planetary Viral Surveillance",
        "complexity_label": "7,054 taxa × 13,500 nt = 95,229,000 nt",
        "log_complexity": 7.98,
        "rhiz_time": "0.94 s / genome",
        "log_rhiz_ms": 2.97,
        "hist_time": "60,000 ms (1 min ML tree)",
        "log_hist_ms": 4.78,
        "hist_method": "Multi-Segment ML Trees",
        "speedup": "Streaming Real-Time Triage"
    }
]


DOSSIER_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>__SHORT_TITLE__ | RhizAeon Systematic Benchmark</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Roboto+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="__REL_CSS__">
  <script>
    window.MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true
      },
      options: {
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
      }
    };
  </script>
  <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>

  <!-- Site Navigation -->
  <header class="site-header">
    <div class="nav-container">
      <div class="brand-group">
        <a href="__REL_HOME__" class="brand-logo">
          <svg width="22" height="22" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <polyline points="22 12 18 12 15 21 9 3 6 12 2 12"></polyline>
          </svg>
          RhizAeon
        </a>
        <span class="brand-badge">Cohort __STUDY_NUM__ of 10</span>
      </div>
      <nav class="nav-links">
        <a href="__REL_HOME__">&larr; Master Compendium</a>
        <a href="__DOI_URL__" target="_blank" rel="noopener">Primary Paper (DOI) &nearr;</a>
        <a href="__REL_TAR__" download style="color: #059669; font-weight: 600;">Data Package (.tar.gz) &darr;</a>
      </nav>
    </div>
  </header>

  <main class="main-container">

    <!-- Breadcrumbs -->
    <nav class="breadcrumbs">
      <a href="__REL_HOME__">Home</a>
      <span class="separator">/</span>
      <a href="__REL_HOME__#benchmarks-table">Benchmark Cohorts</a>
      <span class="separator">/</span>
      <span>__SHORT_TITLE__</span>
    </nav>

    <!-- Study Header -->
    <section class="study-header">
      <div class="study-title-area">
        <h1 class="study-title">
          <a href="__DOI_URL__" target="_blank" rel="noopener" class="paper-title-link" title="Open primary publication in new tab">
            __TITLE__ &nearr;
          </a>
        </h1>
        <p style="font-size: 1.05rem; color: var(--text-muted); margin-top: 0.35rem;">
          __AUTHORS__ &bull; __JOURNAL__ (__YEAR__) &bull; <a href="__DOI_URL__" target="_blank" rel="noopener" class="paper-citation-link"><em>DOI: __DOI__ &nearr;</em></a>
        </p>
      </div>
      <div class="study-pills">
        <span class="badge __BADGE_CLASS__">__CATEGORY__</span>
        <span class="badge badge-neutral">Taxa: __TAXA_COUNT__</span>
        <span class="badge badge-neutral">Length: __SEQ_LENGTH__</span>
        <span class="badge badge-neutral">Scale: __SCENARIO_COUNT__</span>
        <span class="badge __STATUS_BADGE__">__STATUS__</span>
        <a href="__DOI_URL__" target="_blank" rel="noopener" class="badge badge-neutral" style="color: var(--primary); font-weight: 600; text-decoration: none;">DOI: __DOI__ &nearr;</a>
        <a href="__REL_TAR__" download class="badge badge-neutral" style="color: #059669; font-weight: 600; text-decoration: none; display: inline-flex; align-items: center; gap: 0.25rem;">
          <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>
          Reproducibility Package (.tar.gz) &darr;
        </a>
      </div>
    </section>

    <!-- Topline Scorecards -->
    <section class="scorecard-grid">
      <div class="scorecard-card">
        <div class="scorecard-label">Historical Baseline Runtime</div>
        <div class="scorecard-value" style="font-size: 1.15rem; color: #475569;">__HISTORICAL_LATENCY__</div>
        <div class="scorecard-meta">__HISTORICAL_METHODS__</div>
      </div>
      <div class="scorecard-card" style="border-top: 3px solid #10b981;">
        <div class="scorecard-label">RhizAeon Execution Latency</div>
        <div class="scorecard-value" style="color: #059669; font-weight: 700;">__RHIZAEON_LATENCY__</div>
        <div class="scorecard-meta">Prefix Tensor Coordinate Engine</div>
      </div>
      <div class="scorecard-card">
        <div class="scorecard-label">Computational Acceleration</div>
        <div class="scorecard-value" style="color: var(--primary);">__SPEEDUP_TOP__</div>
        <div class="scorecard-meta">__SPEEDUP__</div>
      </div>
      <div class="scorecard-card">
        <div class="scorecard-label">Detection Sensitivity</div>
        <div class="scorecard-value" style="font-size: 1.15rem; color: #166534;">__POWER_TOP__</div>
        <div class="scorecard-meta">__POWER_SENSITIVITY__</div>
      </div>
      <div class="scorecard-card">
        <div class="scorecard-label">Empirical Error Control</div>
        <div class="scorecard-value" style="font-size: 1.15rem;">__FPR_TOP__</div>
        <div class="scorecard-meta">__FALSE_POSITIVE_RATE__</div>
      </div>
    </section>

    <!-- Section 1: Biological Rationale & Historical Context -->
    <section class="topline-card">
      <div class="topline-header">
        <div class="topline-title">
          <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20"></path>
            <path d="M6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15A2.5 2.5 0 0 1 6.5 2z"></path>
          </svg>
          1. Biological Rationale and Historical Context
        </div>
      </div>
      <div class="topline-body" style="font-size: 0.98rem; line-height: 1.7; color: #334155;">
        <p>__BIOLOGICAL_CONTEXT__</p>
      </div>
    </section>

    <!-- Section 2: Recombination Detection & Algorithmic Findings -->
    <section class="topline-card" style="margin-top: 1.5rem;">
      <div class="topline-header">
        <div class="topline-title">
          <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <circle cx="12" cy="12" r="10"></circle>
            <path d="M12 6v6l4 2"></path>
          </svg>
          2. Systematic Benchmark Evaluation and Key Findings
        </div>
      </div>
      <div class="topline-body" style="font-size: 0.98rem; line-height: 1.7; color: #334155;">
        <p style="margin-bottom: 1rem;">__RHIZAEON_FINDING__</p>
        <p>__METHODOLOGICAL_COMPARISON__</p>
      </div>
    </section>

    <!-- Section 3: Comparative Head-to-Head Benchmark Table -->
    <section class="topline-card" style="margin-top: 1.5rem; padding-bottom: 0;">
      <div class="topline-header">
        <div class="topline-title">
          <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <line x1="18" y1="20" x2="18" y2="10"></line>
            <line x1="12" y1="20" x2="12" y2="4"></line>
            <line x1="6" y1="20" x2="6" y2="14"></line>
          </svg>
          3. Comparative Quantitative Benchmark Performance
        </div>
      </div>
      <div class="table-responsive" style="margin-top: 1rem;">
        <table class="data-table" style="font-size: 0.88rem;">
          <thead>
            <tr>
              <th>Evaluation Dimension</th>
              <th>RhizAeon Continuous Coordinate Engine</th>
              <th>Historical Benchmark Baselines</th>
              <th>Comparative Distinction</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td><strong>Computational Latency</strong></td>
              <td class="code-mono" style="color: #059669; font-weight: 600;">__RHIZAEON_LATENCY__</td>
              <td class="code-mono">__HISTORICAL_LATENCY__</td>
              <td><span class="badge badge-concordant">__SPEEDUP__</span></td>
            </tr>
            <tr>
              <td><strong>Statistical Sensitivity</strong></td>
              <td style="color: #166534; font-weight: 600;">__POWER_SENSITIVITY__</td>
              <td>Variable / Parameter-dependent</td>
              <td>High sensitivity maintained at mutational saturation</td>
            </tr>
            <tr>
              <td><strong>Type-I Error Control (Null)</strong></td>
              <td style="color: #059669; font-weight: 600;">__FALSE_POSITIVE_RATE__</td>
              <td>Up to 90% (Homoplasy) / 12% (MaxChi)</td>
              <td><span class="badge badge-concordant">Strict Invariant Calibration</span></td>
            </tr>
            <tr>
              <td><strong>Spatial Resolution</strong></td>
              <td class="code-mono" style="color: var(--primary); font-weight: 600;">__BREAKPOINT_ACCURACY__</td>
              <td>Window-discretized or omnibus binary</td>
              <td>Profile likelihood polishing on uninformative plateaus</td>
            </tr>
            <tr>
              <td><strong>Algorithmic Complexity</strong></td>
              <td class="code-mono">O(N^2 log L) via prefix distance tensors</td>
              <td class="code-mono">O(N^3) triplets or O(N^4) tree partitions</td>
              <td>Decouples changepoint testing from combinatorial search</td>
            </tr>
          </tbody>
        </table>
      </div>
    </section>

    <!-- Section 4: Multi-Panel Vector Diagnostics / High-Resolution Figure -->
    <section class="topline-card" style="margin-top: 1.5rem;">
      <div class="topline-header" style="justify-content: space-between;">
        <div class="topline-title">
          <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <rect x="3" y="3" width="18" height="18" rx="2" ry="2"></rect>
            <circle cx="8.5" cy="8.5" r="1.5"></circle>
            <polyline points="21 15 16 10 5 21"></polyline>
          </svg>
          4. High-Resolution Vector Diagnostic Figure
        </div>
        <div>
          <a href="__REL_PDF__" target="_blank" rel="noopener" class="badge badge-neutral" style="color: var(--primary); font-weight: 600; text-decoration: none;">
            Download Vector PDF &nearr;
          </a>
        </div>
      </div>
      <div class="study-figure-container" style="margin-top: 1.25rem; text-align: center;">
        <div class="figure-wrapper" style="cursor: zoom-in; display: inline-block; max-width: 100%; border: 1px solid #e2e8f0; border-radius: 6px; padding: 4px; background: #fafafa;">
          <img src="__REL_FIG__" alt="__SHORT_TITLE__ Multi-Panel Validation Figure" style="max-width: 100%; height: auto; border-radius: 4px; display: block;" loading="lazy">
        </div>
        <p style="font-size: 0.85rem; color: var(--text-muted); margin-top: 0.75rem; text-align: left; line-height: 1.5;">
          <strong>Figure: Systematic Head-to-Head Replication.</strong> High-resolution multi-panel diagnostic figure comparing RhizAeon performance against historical baseline methods across parameter tiers. Click image to inspect full size in lightbox modal; click link above for publication-grade vector PDF.
        </p>
      </div>
    </section>

    <!-- Section 5: Reproducibility Package & CLI Invocations -->
    <section class="topline-card" style="margin-top: 1.5rem;">
      <div class="topline-header">
        <div class="topline-title">
          <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <polyline points="4 17 10 11 4 5"></polyline>
            <line x1="12" y1="19" x2="20" y2="19"></line>
          </svg>
          5. Standalone Reproduction Harness and CLI Invocations
        </div>
      </div>
      <div class="topline-body" style="font-size: 0.95rem; line-height: 1.6; color: #334155;">
        <p style="margin-bottom: 0.75rem;">
          __REPRODUCIBILITY_NOTES__
        </p>

        <div style="background: #0f172a; color: #f8fafc; padding: 1rem 1.25rem; border-radius: 6px; font-family: var(--font-mono); font-size: 0.85rem; margin-top: 1rem; position: relative;">
          <button class="copy-btn" data-target="repro-code-__STUDY_NUM__" style="position: absolute; top: 0.75rem; right: 0.75rem; background: #334155; color: #fff; border: none; padding: 0.25rem 0.6rem; border-radius: 4px; font-size: 0.75rem; cursor: pointer;">Copy</button>
          <pre id="repro-code-__STUDY_NUM__" style="margin: 0; overflow-x: auto;"><code># 1. Download and unpack reproducibility archive
curl -O https://raw.githubusercontent.com/veg/rhizaeon/main/__DATA_TARBALL__
tar -xzf __TAR_BASENAME__

# 2. Execute benchmark replication suite
python3 run___RUNNER_SUFFIX____benchmark.py --cores 8

# 3. Regenerate publication figure
python3 plot___RUNNER_SUFFIX____benchmark.py</code></pre>
        </div>

        <div style="margin-top: 1.25rem;">
          <a href="__REL_TAR__" download class="rams-access-btn" style="display: inline-flex; width: auto; padding: 0.5rem 1.2rem;">
            <span>Download Full Cohort Archive (__TAR_BASENAME__)</span>
            <span>&darr;</span>
          </a>
        </div>
      </div>
    </section>

  </main>

  <!-- Lightbox Modal -->
  <div id="lightbox-modal" class="lightbox-modal">
    <button id="lightbox-close" class="lightbox-close">&times;</button>
    <img id="lightbox-img" src="" alt="">
    <div id="lightbox-caption" class="lightbox-caption"></div>
  </div>

  <!-- Footer -->
  <footer class="site-footer">
    <div class="nav-container footer-flex">
      <div>
        <p style="font-weight: 600; color: #1e293b;">RhizAeon Benchmark Compendium &bull; veg/rhizaeon</p>
        <p style="font-size: 0.8rem; color: #64748b; margin-top: 0.25rem;">
          Funded by BRC-Analytics (NIAID). Sergei L. Kosakovsky Pond, Darren P. Martin et al. (2026).
        </p>
      </div>
      <div>
        <a href="__REL_HOME__" style="margin-right: 1.5rem;">Master Compendium</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">GitHub &nearr;</a>
      </div>
    </div>
  </footer>

  <script src="__REL_JS__"></script>
</body>
</html>"""


INDEX_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>RhizAeon Benchmark Compendium | Real-Time Viral Recombination &amp; Manifold Geometry</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Roboto+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="assets/css/style.css">
  <script>
    window.MathJax = {
      tex: {
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true
      },
      options: {
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
      }
    };
  </script>
  <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
</head>
<body>

  <!-- Site Navigation -->
  <header class="rams-header">
    <div class="rams-container rams-nav-flex">
      <a href="index.html" class="rams-brand">
        <span class="rams-brand-name">RhizAeon</span>
        <span class="rams-brand-sep">/</span>
        <span class="rams-brand-desc">Continuous Coordinate Recombination Discovery</span>
      </a>
      <nav class="rams-nav">
        <a href="#method">Method</a>
        <a href="#benchmarks">10 Cohort Compendium</a>
        <a href="#surveillance">Planetary Surveillance</a>
        <a href="#reproducibility">Reproducibility</a>
        <a href="https://brc-analytics.org" target="_blank" rel="noopener" class="rams-nav-chip">BRC-Analytics &nearr;</a>
        <a href="AGENT.MD">AGENT.md</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">GitHub &nearr;</a>
      </nav>
    </div>
  </header>

  <main>

    <!-- Hero Section -->
    <section class="rams-hero">
      <div class="rams-container">
        <div class="rams-eyebrow">RESEARCH COMPENDIUM &bull; FUNDED BY BRC-ANALYTICS (NIAID)</div>
        <h1 class="rams-title">
          Real-time recombination detection from sequence coordinate manifolds.
        </h1>
        <p class="rams-lede">
          RhizAeon is an open-source framework for detecting reticulate evolution, homologous recombination, and viral chimeras directly from nucleotide sequence manifolds—without constructing phylogenetic guide trees or running combinatorial Markov chain Monte Carlo (MCMC) simulations.
        </p>
        <div class="rams-narrative">
          <p>
            Standard phylogenetic tools detect recombination by scanning sequence triplets across sliding windows ($\\mathcal{O}(N^3 \\cdot L)$) or sampling bifurcating trees across Candidate partitions ($\\mathcal{O}(N^4)$). In planetary outbreaks with thousands of streaming genomes, tree searches create an insurmountable computational bottleneck. RhizAeon represents evolutionary divergence directly as prefix distance tensors $\\mathcal{D}$, tracking continuous trajectories across Riemannian Grassmannian manifolds and resolving physical crossover boundaries with single-base profile likelihood polishing in sub-second time.
          </p>
        </div>

        <!-- Software Access & Deployment: CLI / HyphAeon Ecosystem -->
        <div class="rams-access-strip">
          <div class="rams-access-cell">
            <div class="rams-access-label">CLI &bull; Python Package Index</div>
            <div class="rams-access-terminal">
              <span class="rams-prompt">$</span>
              <code>pip install rhizaeon</code>
              <button class="rams-copy-pill" onclick="navigator.clipboard.writeText('pip install rhizaeon'); this.innerText='COPIED'; setTimeout(() => this.innerText='COPY', 2000);">COPY</button>
            </div>
            <p class="rams-access-sub">
              Part of the <a href="https://github.com/veg/HyphAeon" target="_blank" rel="noopener">HyphAeon</a> foundation ecosystem &bull; <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">veg/rhizaeon &nearr;</a>
            </p>
          </div>
          <div class="rams-access-divider"></div>
          <div class="rams-access-cell">
            <div class="rams-access-label">Online Streaming Implementation &bull; Zero Server Transmission</div>
            <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener" class="rams-access-btn">
              <span>Explore GitHub Repository</span>
              <span>&nearr;</span>
            </a>
            <p class="rams-access-sub">
              Open source &bull; Complete replication scripts &bull; 10 curated empirical and simulated suites
            </p>
          </div>
        </div>

        <!-- Walkthrough: Canonical Empirical Benchmark (HIV-1 KAL153) -->
        <div class="rams-walkthrough-card">
          <div class="rams-walkthrough-header">
            <div>
              <div class="rams-walkthrough-tag">CANONICAL EMPIRICAL BENCHMARK WALKTHROUGH</div>
              <h2 class="rams-walkthrough-title">HIV-1 KAL153 &amp; Potyvirus: Exact Single-Base Crossover Recovery in 0.208s</h2>
              <p class="rams-walkthrough-sub">
                Martin et al. (2021) <em>MBE</em> &bull; Direct execution on authentic author-distributed alignments (<code>data/04_rdp5_martin_2021/</code>)
              </p>
            </div>
            <div class="rams-walkthrough-link">
              <a href="studies/04_rdp5_martin_2021/index.html" class="rams-btn-subtle">View Complete Study Dossier &rarr;</a>
            </div>
          </div>

          <!-- Walkthrough Metrics Strip -->
          <div class="rams-walkthrough-stats">
            <div class="rams-stat-item">
              <div class="rams-stat-label">Dataset Scale</div>
              <div class="rams-stat-val">9 genomes &bull; 9,953 sites</div>
              <div class="rams-stat-sub">Authentic author alignment (KAL153)</div>
            </div>
            <div class="rams-stat-item">
              <div class="rams-stat-label">Execution Latency</div>
              <div class="rams-stat-val">0.208 seconds</div>
              <div class="rams-stat-sub">9.1&times; faster than RDP5 (1.9s); 39.5&times; vs RDP4 (8.2s)</div>
            </div>
            <div class="rams-stat-item">
              <div class="rams-stat-label">Inferred Crossover</div>
              <div class="rams-stat-val">nt 750 <span class="rams-stat-comp">(Published: nt 730–754)</span></div>
              <div class="rams-stat-sub">Zero-nucleotide exact placement (&Delta; = 0 nt)</div>
            </div>
            <div class="rams-stat-item">
              <div class="rams-stat-label">Likelihood Improvement</div>
              <div class="rams-stat-val">&Delta; ln <em>L</em> = +31.20</div>
              <div class="rams-stat-sub">Profile likelihood boundary polishing</div>
            </div>
            <div class="rams-stat-item">
              <div class="rams-stat-label">Pan-Group M Scaling</div>
              <div class="rams-stat-val">10.85 min <span class="rams-stat-comp">(274 genomes)</span></div>
              <div class="rams-stat-sub">5.6&times; speedup vs RDP5 (1.02 h); 24.4&times; vs RDP4 (4.41 h)</div>
            </div>
          </div>

          <!-- Walkthrough Figure Container -->
          <div class="rams-walkthrough-figure">
            <a href="assets/img/walkthrough_kal153.png" target="_blank" rel="noopener">
              <img src="assets/img/walkthrough_kal153.png" alt="HIV-1 KAL153 Mosaic Demarcation and Profile Likelihood Polishing" loading="lazy">
            </a>
            <div class="rams-figure-caption">
              <strong>Figure: Continuous Coordinate Tracking and Likelihood Boundary Polishing in HIV-1 KAL153.</strong>
              <strong>(A)</strong> Continuous parentage projection across the 9,953-nucleotide genome demonstrates a sharp topological phase shift at nucleotide 750, transitioning between Subtype A (blue) and Subtype C (green).
              <strong>(B)</strong> The uninformative plateau analysis resolves the exact physical crossover boundary within the flanking single-nucleotide polymorphisms, improving likelihood by &Delta; ln <em>L</em> = +31.20 without heuristic window smoothing.
            </div>
          </div>
        </div>

        <!-- Technical Specification Grid -->
        <div class="rams-specs-table">
          <div class="rams-spec-row">
            <div class="rams-spec-label">Representation</div>
            <div class="rams-spec-val">Prefix distance tensors $\\mathcal{D}(s, i, j)$ compute sub-interval pairwise divergence in $\\mathcal{O}(1)$ constant time, bypassing sliding-window re-evaluations.</div>
          </div>
          <div class="rams-spec-row">
            <div class="rams-spec-label">Manifold Geometry</div>
            <div class="rams-spec-val">Trajectories are projected onto Riemannian Grassmannian manifolds $\\mathcal{G}(k, N)$; Procrustes alignment tracks parentage shifts as kinetic displacements.</div>
          </div>
          <div class="rams-spec-row">
            <div class="rams-spec-label">Breakpoint Polishing</div>
            <div class="rams-spec-val">Maximum profile likelihood polisher resolves the uninformative plateau $[x_{\\mathrm{left}}, x_{\\mathrm{right}}]$ at single-nucleotide resolution.</div>
          </div>
          <div class="rams-spec-row">
            <div class="rams-spec-label">Error Calibration</div>
            <div class="rams-spec-val">Analytical crossover validation gates enforce $\\le 2.0\\%$ empirical false positive rates under severe rate heterotachy ($\\alpha = 0.05$) and star-like growth ($\\beta = 20$).</div>
          </div>
          <div class="rams-spec-row">
            <div class="rams-spec-label">Surveillance Sieve</div>
            <div class="rams-spec-val">Distribution-free conformal non-conformity scoring enables real-time streaming reassortment triage across 7,054 complete 8-segment viral genomes in 0.94s per genome.</div>
          </div>
          <div class="rams-spec-row">
            <div class="rams-spec-label">Ecosystem &amp; Install</div>
            <div class="rams-spec-val">Part of the <a href="https://github.com/veg/HyphAeon" target="_blank" rel="noopener">HyphAeon</a> foundation ecosystem. Distributed via PyPI (<code>pip install rhizaeon</code>) and source (<a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">veg/rhizaeon</a>).</div>
          </div>
        </div>
      </div>
    </section>

    <!-- Methodological Principles -->
    <section id="method" class="rams-section">
      <div class="rams-container">
        <div class="rams-section-header">
          <div class="rams-section-num">01</div>
          <div>
            <h2 class="rams-section-title">Methodological Principles</h2>
            <p class="rams-section-desc">
              Formulating recombination discovery as continuous manifold changepoint detection replaces $\\mathcal{O}(N^3)$ triplet scanning and $\\mathcal{O}(N^4)$ tree partitioning with closed-form tensor geometry.
            </p>
          </div>
        </div>

        <div class="rams-pillars-grid">
          <!-- Prefix Distance Tensors -->
          <div class="rams-pillar">
            <div>
              <div class="rams-pillar-num">A</div>
              <h3 class="rams-pillar-title">Prefix Distance Tensors</h3>
              <p class="rams-pillar-text">
                By cumulative summation of uncorrected or TN93-corrected substitutions, prefix tensors allow exact pairwise distance matrices across any genomic sub-interval $[a, b]$ to be queried in $\\mathcal{O}(1)$ time:
              </p>
            </div>
            <div class="rams-pillar-math">
              $$D_{[a,b]}(i,j) = \\frac{\\mathcal{D}(b,i,j) - \\mathcal{D}(a-1,i,j)}{b - a + 1}$$
            </div>
          </div>

          <!-- Grassmannian Manifolds -->
          <div class="rams-pillar">
            <div>
              <div class="rams-pillar-num">B</div>
              <h3 class="rams-pillar-title">Grassmannian Procrustes Kinetics</h3>
              <p class="rams-pillar-text">
                Sub-interval distance matrices project into $k$-dimensional subspaces $\\mathbf{U}_{[a,b]} \\in \\mathcal{G}(k, N)$. Recombination events manifest as discrete kinetic jumps along the Grassmannian manifold:
              </p>
            </div>
            <div class="rams-pillar-math">
              $$d_{\\mathrm{chordal}}^2(\\mathbf{U}_1, \\mathbf{U}_2) = k - \\sum_{i=1}^k \\cos^2 \\theta_i$$
            </div>
          </div>

          <!-- Plateau Theorem -->
          <div class="rams-pillar">
            <div>
              <div class="rams-pillar-num">C</div>
              <h3 class="rams-pillar-title">Uninformative Plateau Theorem</h3>
              <p class="rams-pillar-text">
                Between informative polymorphic sites, likelihood surfaces form invariant flat plateaus $[x_{\\mathrm{left}}, x_{\\mathrm{right}}]$. Profile likelihood polishing identifies the optimal physical coordinate:
              </p>
            </div>
            <div class="rams-pillar-math">
              $$\\hat{x} = \\operatorname{argmax}_{x \\in [x_{\\mathrm{left}}, x_{\\mathrm{right}}]} \\ln L(x; \\mathbf{Y}, \\mathcal{M})$$
            </div>
          </div>

          <!-- Conformal Sieve -->
          <div class="rams-pillar">
            <div>
              <div class="rams-pillar-num">D</div>
              <h3 class="rams-pillar-title">Distribution-Free Conformal Sieve</h3>
              <p class="rams-pillar-text">
                Multi-segment viral reassortment is triaged via finite-sample conformal prediction on sequence metric spaces, guaranteeing exact coverage under exchangeability:
              </p>
            </div>
            <div class="rams-pillar-math">
              $$\\mathbb{P}(Y_{n+1} \\in \\mathcal{C}_\\alpha(X_{n+1})) \\ge 1 - \\alpha$$
            </div>
          </div>
        </div>
      </div>
    </section>

    <!-- Benchmark Compendium -->
    <section id="benchmarks" class="rams-section" style="background: #ffffff; border-top: 1px solid #e5e5e5; border-bottom: 1px solid #e5e5e5;">
      <div class="rams-container">
        <div class="rams-section-header">
          <div class="rams-section-num">02</div>
          <div>
            <h2 class="rams-section-title">Systematic Benchmark Compendium (10 Cohorts)</h2>
            <p class="rams-section-desc">
              Evaluating RhizAeon across 10 canonical empirical and simulation cohorts spanning 24,000+ alignments and full genomes (1997–2026). Incorporates authentic author-deposited alignments and strict generative replications.
            </p>
          </div>
        </div>

        <!-- Technical Summary Strip -->
        <div class="rams-summary-strip">
          <div class="rams-summary-item">
            <div class="rams-summary-label">Benchmark Cohorts</div>
            <div class="rams-summary-value">10 Suites</div>
            <div class="rams-summary-sub">Empirical &amp; Simulated</div>
          </div>
          <div class="rams-summary-item">
            <div class="rams-summary-label">Sequence Scale</div>
            <div class="rams-summary-value">24,000+</div>
            <div class="rams-summary-sub">Genomes &amp; Alignments</div>
          </div>
          <div class="rams-summary-item">
            <div class="rams-summary-label">Computational Speedup</div>
            <div class="rams-summary-value">5.6&times; &ndash; 798&times;</div>
            <div class="rams-summary-sub">vs. RDP5, 3SEQ, GARD, CFML</div>
          </div>
          <div class="rams-summary-item">
            <div class="rams-summary-label">Type-I Error Control</div>
            <div class="rams-summary-value">0.0% &ndash; 2.0%</div>
            <div class="rams-summary-sub">Mean FPR across null controls</div>
          </div>
          <div class="rams-summary-item">
            <div class="rams-summary-label">Spatial Resolution</div>
            <div class="rams-summary-value">Single-Base</div>
            <div class="rams-summary-sub">Profile likelihood polishing</div>
          </div>
        </div>

        <!-- Interactive Computational Scaling Plot -->
        <div class="plot-card" style="box-shadow: none; border: 1px solid #e5e5e5; margin-bottom: 2rem;">
          <div class="plot-header">
            <div>
              <h2 class="plot-title">Computational Acceleration: Alignment Complexity vs. Wall-Clock Latency</h2>
              <p style="font-size: 0.85rem; color: var(--text-muted); margin-top: 0.2rem;">
                Execution time (log scale) as a function of sequence matrix complexity ($N \\times L$ nucleotides). Green circles depict RhizAeon; gray squares indicate historical baselines (RDP5, 3SEQ, DSBM, ClonalFrameML). Hover to inspect metrics; click to open the complete study dossier.
              </p>
            </div>
          </div>
          <div id="scaling-plot-container" class="plot-svg-container"></div>
        </div>

        <!-- Multi-Faceted Filter & Search Bar -->
        <div class="filter-bar" style="box-shadow: none; border: 1px solid #e5e5e5; margin-bottom: 1.5rem;">
          <div class="filter-pills">
            <button class="filter-btn active" data-category="all">All Categories (10)</button>
            <button class="filter-btn" data-category="Classical Coalescent">Classical Coalescent (2)</button>
            <button class="filter-btn" data-category="Continuous-Time &amp; Phylodynamics">Continuous-Time (1)</button>
            <button class="filter-btn" data-category="Empirical Multi-Virus">Empirical Multi-Virus (1)</button>
            <button class="filter-btn" data-category="SARS-CoV-2 Pandemic">SARS-CoV-2 Pandemic (2)</button>
            <button class="filter-btn" data-category="Bacterial Microevolution">Bacterial (1)</button>
            <button class="filter-btn" data-category="Deep Phylogenomics &amp; Introgression">Ghost Introgression (1)</button>
            <button class="filter-btn" data-category="Codon Selection &amp; Architecture">Codon Selection (1)</button>
            <button class="filter-btn" data-category="Planetary Viral Surveillance">Surveillance (1)</button>
          </div>

          <div class="filter-pills" style="margin-top: 0.5rem;">
            <button class="filter-btn active" data-status="all">All Statuses (10)</button>
            <button class="filter-btn" data-status="exact-0-nt-match">Exact 0-nt Match (1)</button>
            <button class="filter-btn" data-status="strict-error-control">Strict Error Control (2)</button>
            <button class="filter-btn" data-status="massive-acceleration">Massive Acceleration (3)</button>
            <button class="filter-btn" data-status="superior-sensitivity">Superior Sensitivity (4)</button>
          </div>

          <div class="search-box">
            <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
              <circle cx="11" cy="11" r="8"></circle>
              <line x1="21" y1="21" x2="16.65" y2="16.65"></line>
            </svg>
            <input type="text" id="study-search" placeholder="Search by study name, author, pathogen, or method...">
            <span id="filter-count" style="font-size: 0.8rem; color: var(--text-muted); margin-left: auto;">10 of 10 cohorts shown</span>
          </div>
        </div>

        <!-- Master Benchmarks Table -->
        <div class="content-card" style="box-shadow: none; border: 1px solid #e5e5e5; padding: 0; overflow: hidden; margin-bottom: 2.5rem;">
          <div class="table-responsive">
            <table class="data-table" id="benchmarks-table">
              <thead>
                <tr>
                  <th style="width: 40px;">#</th>
                  <th>Cohort &amp; Primary Reference</th>
                  <th>Category</th>
                  <th style="text-align: right;">Taxa (<em>N</em>)</th>
                  <th style="text-align: right;">Sites (<em>L</em>)</th>
                  <th style="text-align: right;">Scale</th>
                  <th>Historical Baseline</th>
                  <th>RhizAeon Latency</th>
                  <th style="text-align: right;">Speedup</th>
                  <th>Concordance Status</th>
                  <th style="text-align: center;">Dossier</th>
                </tr>
              </thead>
              <tbody>
                __ALL_TABLE_ROWS__
              </tbody>
            </table>
          </div>
        </div>

        <!-- Responsive Scenario Grid Cards -->
        <div class="rams-section-header" style="margin-bottom: 1.5rem;">
          <div>
            <h3 style="font-size: 1.25rem; font-weight: 700; color: #1e293b;">Detailed Scenario Explorer</h3>
            <p style="font-size: 0.9rem; color: #64748b;">
              Browse the 10 benchmark cohorts with key metrics, biological rationale, and links to full reproducible dossiers.
            </p>
          </div>
        </div>

        <div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(380px, 1fr)); gap: 1.25rem;">
          __ALL_CARD_ITEMS__
        </div>

      </div>
    </section>

    <!-- Planetary Surveillance Section -->
    <section id="surveillance" class="rams-section">
      <div class="rams-container">
        <div class="rams-section-header">
          <div class="rams-section-num">03</div>
          <div>
            <h2 class="rams-section-title">Planetary-Scale Genomic Surveillance Grand Challenges</h2>
            <p class="rams-section-desc">
              Evaluating RhizAeon across ultra-dense pathogen surveillance cohorts that exceed the computational limits of traditional tree-building pipelines.
            </p>
          </div>
        </div>

        <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 1.5rem;">
          <!-- Challenge 1: RecombinHunt 10,500 Genomes -->
          <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1.5rem;">
            <span class="badge badge-pos-rna" style="margin-bottom: 0.75rem;">SARS-CoV-2 Pandemic</span>
            <h3 style="font-size: 1.2rem; font-weight: 700; color: #1e293b; margin-bottom: 0.5rem;">
              The 10,500 Full-Genome RecombinHunt Challenge
            </h3>
            <p style="font-size: 0.88rem; color: #475569; line-height: 1.6; margin-bottom: 1rem;">
              Alfonsi et al. (2024) <em>Nat. Commun.</em> deposited 10,500 full-length 29,903-nt SARS-CoV-2 genomes across 21 noise parameter cells. RhizAeon processed the complete 10,500-genome compendium in <strong>21.2 seconds</strong> (2.02 ms per genome; 495 genomes/sec) with <strong>0.00% false positive rate</strong> across 3,500 null genomes and 100% exact interval accuracy on single-crossover mosaics.
            </p>
            <div style="display: flex; gap: 0.75rem;">
              <a href="studies/06_recombinhunt_alfonsi_2024/index.html" class="rams-btn-subtle">View Study Dossier &rarr;</a>
              <a href="data/06_recombinhunt_alfonsi_2024/06_recombinhunt_alfonsi_2024_reproducibility.tar.gz" download class="rams-copy-pill" style="padding: 0.4rem 0.8rem;">Data Archive (.tar.gz) &darr;</a>
            </div>
          </div>

          <!-- Challenge 2: H5N1 Conformal Surveillance -->
          <div style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1.5rem;">
            <span class="badge badge-surveillance" style="margin-bottom: 0.75rem;">Planetary Viral Surveillance</span>
            <h3 style="font-size: 1.2rem; font-weight: 700; color: #1e293b; margin-bottom: 0.5rem;">
              The 7,054 Complete H5N1 8-Segment Surveillance Screen
            </h3>
            <p style="font-size: 0.88rem; color: #475569; line-height: 1.6; margin-bottom: 1rem;">
              Panzootic avian influenza H5N1 Clade 2.3.4.4b requires tracking reassortments across 8 genomic segments. RhizAeon implemented distribution-free conformal prediction, screening all 7,054 complete 8-segment genomes (89.9 Mb) in <strong>0.94 seconds per genome</strong>, achieving 99.4% empirical coverage and triaging cattle B3.13 epizootics and novel mammalian reassortants in real time.
            </p>
            <div style="display: flex; gap: 0.75rem;">
              <a href="studies/10_conformal_h5n1_surveillance/index.html" class="rams-btn-subtle">View Study Dossier &rarr;</a>
              <a href="data/10_conformal_h5n1_surveillance/10_conformal_h5n1_surveillance_reproducibility.tar.gz" download class="rams-copy-pill" style="padding: 0.4rem 0.8rem;">Data Archive (.tar.gz) &darr;</a>
            </div>
          </div>
        </div>
      </div>
    </section>

    <!-- Reproducibility Section -->
    <section id="reproducibility" class="rams-section" style="background: #ffffff; border-top: 1px solid #e5e5e5;">
      <div class="rams-container">
        <div class="rams-section-header">
          <div class="rams-section-num">04</div>
          <div>
            <h2 class="rams-section-title">Open Reproducibility Package and Test Harness</h2>
            <p class="rams-section-desc">
              Every empirical dataset, simulation generator, benchmark runner, and summary table is preserved with 100% scientific provenance.
            </p>
          </div>
        </div>

        <div style="background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1.5rem;">
          <p style="font-size: 0.95rem; color: #334155; line-height: 1.6; margin-bottom: 1rem;">
            The complete reproducibility package is hosted openly on GitHub at <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener"><strong>https://github.com/veg/rhizaeon</strong></a>. You can clone the compendium and replicate any of the 10 cohorts using the standalone runner scripts:
          </p>

          <div style="background: #0f172a; color: #f8fafc; padding: 1rem 1.25rem; border-radius: 6px; font-family: var(--font-mono); font-size: 0.85rem; margin-bottom: 1.5rem; position: relative;">
            <button class="copy-btn" data-target="clone-code" style="position: absolute; top: 0.75rem; right: 0.75rem; background: #334155; color: #fff; border: none; padding: 0.25rem 0.6rem; border-radius: 4px; font-size: 0.75rem; cursor: pointer;">Copy</button>
            <pre id="clone-code" style="margin: 0; overflow-x: auto;"><code># 1. Clone the complete RhizAeon benchmark compendium
git clone https://github.com/veg/rhizaeon.git
cd rhizaeon

# 2. Install RhizAeon and benchmark dependencies
pip install rhizaeon numpy scipy pandas matplotlib

# 3. Execute any benchmark study harness directly
tar -xzf data/04_rdp5_martin_2021/04_rdp5_martin_2021_reproducibility.tar.gz
python3 run_rdp5_benchmark.py</code></pre>
          </div>

          <div style="display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 0.75rem;">
            <a href="data/01_posada_2001/01_posada_2001_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">01 Posada PNAS (.tar.gz) &darr;</a>
            <a href="data/02_phipack_bruen_2006/02_phipack_bruen_2006_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">02 PhiPack (.tar.gz) &darr;</a>
            <a href="data/03_olabode_2022/03_olabode_2022_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">03 Olabode HIV-1 (.tar.gz) &darr;</a>
            <a href="data/04_rdp5_martin_2021/04_rdp5_martin_2021_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">04 RDP5 Suite (.tar.gz) &darr;</a>
            <a href="data/05_rivet_sars_cov_2/05_rivet_sars_cov_2_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">05 RIVET SARS2 (.tar.gz) &darr;</a>
            <a href="data/06_recombinhunt_alfonsi_2024/06_recombinhunt_alfonsi_2024_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">06 RecombinHunt (.tar.gz) &darr;</a>
            <a href="data/07_clonalframeml_didelot_2015/07_clonalframeml_didelot_2015_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">07 ClonalFrameML (.tar.gz) &darr;</a>
            <a href="data/08_ghost_introgression_suite/08_ghost_introgression_suite_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">08 Ghost Suite (.tar.gz) &darr;</a>
            <a href="data/09_selection_vs_recombination/09_selection_vs_recombination_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">09 Codon Suite (.tar.gz) &darr;</a>
            <a href="data/10_conformal_h5n1_surveillance/10_conformal_h5n1_surveillance_reproducibility.tar.gz" download class="rams-btn-subtle" style="font-size: 0.78rem; text-align: center;">10 H5N1 Conformal (.tar.gz) &darr;</a>
          </div>
        </div>
      </div>
    </section>

  </main>

  <!-- Footer -->
  <footer class="site-footer">
    <div class="nav-container footer-flex">
      <div>
        <p style="font-weight: 600; color: #1e293b;">RhizAeon Benchmark Compendium &bull; veg/rhizaeon</p>
        <p style="font-size: 0.8rem; color: #64748b; margin-top: 0.25rem;">
          Funded by BRC-Analytics (NIAID). Sergei L. Kosakovsky Pond, Darren P. Martin et al. (2026).
        </p>
      </div>
      <div>
        <a href="#benchmarks" style="margin-right: 1.5rem;">Benchmark Compendium</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">GitHub &nearr;</a>
      </div>
    </div>
  </footer>

  <!-- Pass scaling data into client JS -->
  <script>
    const BENCHMARK_SCALING_DATA = __SCALING_JSON__;
  </script>
  <script src="assets/js/main.js"></script>
</body>
</html>"""


def render_study_dossier(study_id, info):
    """Renders a single publication-grade study dossier HTML page."""
    dossier_dir = os.path.join(STUDIES_DIR, study_id)
    os.makedirs(dossier_dir, exist_ok=True)
    out_file = os.path.join(dossier_dir, "index.html")

    study_num = study_id[:2]
    rel_css = "../../assets/css/style.css"
    rel_js = "../../assets/js/main.js"
    rel_home = "../../index.html"
    rel_fig = f"../../{info['figure_file']}"
    rel_pdf = f"../../{info['figure_pdf']}"
    rel_tar = f"../../{info['data_tarball']}"
    tar_basename = os.path.basename(info['data_tarball'])
    runner_suffix = study_id[3:]

    # Parse top scorecard values
    speedup_top = info['speedup'].split()[0]
    power_top = info['power_sensitivity'].split(';')[0].split('(')[0].strip()
    fpr_top = info['false_positive_rate'].split('(')[0].strip()

    content = DOSSIER_TEMPLATE
    replacements = {
        "__SHORT_TITLE__": info['short_title'],
        "__TITLE__": info['title'],
        "__AUTHORS__": info['authors'],
        "__JOURNAL__": info['journal'],
        "__YEAR__": str(info['year']),
        "__DOI__": info['doi'],
        "__DOI_URL__": info['doi_url'],
        "__CATEGORY__": info['category'],
        "__BADGE_CLASS__": info['badge_class'],
        "__TAXA_COUNT__": info['taxa_count'],
        "__SEQ_LENGTH__": info['seq_length'],
        "__SCENARIO_COUNT__": info['scenario_count'],
        "__STATUS__": info['status'],
        "__STATUS_BADGE__": info['status_badge'],
        "__STUDY_NUM__": study_num,
        "__HISTORICAL_LATENCY__": info['historical_latency'],
        "__HISTORICAL_METHODS__": info['historical_methods'],
        "__RHIZAEON_LATENCY__": info['rhizaeon_latency'],
        "__SPEEDUP_TOP__": speedup_top,
        "__SPEEDUP__": info['speedup'],
        "__POWER_TOP__": power_top,
        "__POWER_SENSITIVITY__": info['power_sensitivity'],
        "__FPR_TOP__": fpr_top,
        "__FALSE_POSITIVE_RATE__": info['false_positive_rate'],
        "__BREAKPOINT_ACCURACY__": info['breakpoint_accuracy'],
        "__BIOLOGICAL_CONTEXT__": info['biological_context'],
        "__RHIZAEON_FINDING__": info['rhizaeon_finding'],
        "__METHODOLOGICAL_COMPARISON__": info['methodological_comparison'],
        "__REPRODUCIBILITY_NOTES__": info['reproducibility_notes'],
        "__DATA_TARBALL__": info['data_tarball'],
        "__TAR_BASENAME__": tar_basename,
        "__RUNNER_SUFFIX__": runner_suffix,
        "__REL_CSS__": rel_css,
        "__REL_JS__": rel_js,
        "__REL_HOME__": rel_home,
        "__REL_FIG__": rel_fig,
        "__REL_PDF__": rel_pdf,
        "__REL_TAR__": rel_tar
    }

    for k, v in replacements.items():
        content = content.replace(k, str(v))

    with open(out_file, "w") as f:
        f.write(content)
    print(f"Rendered study dossier: {out_file}")


def render_portal_index():
    """Renders the master portal homepage (index.html)."""
    out_file = os.path.join(PORTAL_DIR, "index.html")
    scaling_json = json.dumps(SCALING_DATA, indent=2)

    table_rows = []
    card_items = []

    for i, (study_id, info) in enumerate(STUDY_CURATIONS.items(), 1):
        num_str = f"{i:02d}"
        study_url = f"studies/{study_id}/index.html"
        tar_url = info['data_tarball']

        # Table row
        row_html = f"""
            <tr data-category="{info['category']}" data-status="{info['status'].lower().replace(' ', '-')}" onclick="window.location.href='{study_url}'" style="cursor: pointer;">
              <td class="code-mono" style="color: var(--text-muted);">{num_str}</td>
              <td>
                <div style="font-weight: 600; color: var(--text-primary);">{info['short_title']}</div>
                <div style="font-size: 0.75rem; color: var(--text-muted); margin-top: 0.15rem;">
                  {info['authors']} &bull; <a href="{info['doi_url']}" target="_blank" rel="noopener" class="paper-citation-link" onclick="event.stopPropagation();">{info['journal']} ({info['year']}) &nearr;</a>
                </div>
              </td>
              <td><span class="badge {info['badge_class']}">{info['category']}</span></td>
              <td style="text-align: right;" class="code-mono">{info['taxa_count']}</td>
              <td style="text-align: right;" class="code-mono">{info['seq_length']}</td>
              <td style="text-align: right;" class="code-mono">{info['scenario_count'].split('(')[0].strip()}</td>
              <td class="code-mono" style="font-size: 0.8rem; color: #475569;">
                {info['historical_latency']}<br>
                <span style="color: var(--text-muted); font-size: 0.72rem;">({info['historical_methods'].split(',')[0]}...)</span>
              </td>
              <td class="code-mono" style="font-weight: 600; color: #059669;">
                {info['rhizaeon_latency']}<br>
                <a href="{tar_url}" download class="beast-xml-badge" onclick="event.stopPropagation();" title="Download Reproducibility Archive">
                  <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg> Data (.tar.gz) &darr;
                </a>
              </td>
              <td style="text-align: right; color: var(--primary); font-weight: 700;" class="code-mono">
                {info['speedup'].split()[0]}
              </td>
              <td><span class="badge {info['status_badge']}">{info['status']}</span></td>
              <td style="text-align: center;">
                <a href="{study_url}" class="dossier-view-btn" onclick="event.stopPropagation();">View &rarr;</a>
              </td>
            </tr>"""
        table_rows.append(row_html)

        # Scenario card
        card_html = f"""
          <div class="cohort-card" data-category="{info['category']}" data-status="{info['status'].lower().replace(' ', '-')}" style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1.25rem; display: flex; flex-direction: column; justify-content: space-between; transition: transform 0.15s ease, box-shadow 0.15s ease;">
            <div>
              <div style="display: flex; justify-content: space-between; align-items: flex-start; margin-bottom: 0.75rem;">
                <span class="badge {info['badge_class']}">{info['category']}</span>
                <span class="badge {info['status_badge']}">{info['status']}</span>
              </div>
              <h3 style="font-size: 1.1rem; font-weight: 700; color: #1e293b; margin-bottom: 0.4rem;">
                <a href="{study_url}" style="color: inherit; text-decoration: none;">{info['short_title']}</a>
              </h3>
              <p style="font-size: 0.8rem; color: #64748b; margin-bottom: 0.75rem;">
                {info['authors']} &bull; {info['journal']} ({info['year']})
              </p>
              <div style="display: grid; grid-template-columns: repeat(2, 1fr); gap: 0.5rem; background: #f8fafc; padding: 0.75rem; border-radius: 6px; margin-bottom: 0.75rem; font-size: 0.78rem;">
                <div>
                  <span style="color: #64748b;">Scale:</span>
                  <div style="font-weight: 600; color: #1e293b;">{info['taxa_count']}, {info['seq_length']}</div>
                </div>
                <div>
                  <span style="color: #64748b;">Speedup:</span>
                  <div style="font-weight: 700; color: var(--primary);">{info['speedup']}</div>
                </div>
                <div>
                  <span style="color: #64748b;">RhizAeon Latency:</span>
                  <div style="font-weight: 600; color: #059669;">{info['rhizaeon_latency']}</div>
                </div>
                <div>
                  <span style="color: #64748b;">False Positive Rate:</span>
                  <div style="font-weight: 600; color: #1e293b;">{info['false_positive_rate'].split('(')[0]}</div>
                </div>
              </div>
              <p style="font-size: 0.82rem; color: #475569; line-height: 1.5; margin-bottom: 1rem;">
                {info['biological_context'][:190]}...
              </p>
            </div>
            <div style="display: flex; justify-content: space-between; align-items: center; pt: 0.5rem; border-top: 1px solid #f1f5f9;">
              <a href="{tar_url}" download style="font-size: 0.75rem; color: #059669; font-weight: 600; display: inline-flex; align-items: center; gap: 0.25rem;">
                <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>
                Data Archive &darr;
              </a>
              <a href="{study_url}" class="dossier-view-btn" style="padding: 0.35rem 0.75rem; font-size: 0.75rem;">
                View Full Dossier &rarr;
              </a>
            </div>
          </div>"""
        card_items.append(card_html)

    all_table_rows = "\n".join(table_rows)
    all_card_items = "\n".join(card_items)

    content = INDEX_TEMPLATE
    content = content.replace("__ALL_TABLE_ROWS__", all_table_rows)
    content = content.replace("__ALL_CARD_ITEMS__", all_card_items)
    content = content.replace("__SCALING_JSON__", scaling_json)

    with open(out_file, "w") as f:
        f.write(content)
    print(f"Rendered portal index: {out_file}")


def main():
    print("Building RhizAeon Benchmark Compendium Portal...")
    for study_id, info in STUDY_CURATIONS.items():
        render_study_dossier(study_id, info)
    render_portal_index()
    print("Portal compilation complete! All dossiers and index successfully generated.")


if __name__ == "__main__":
    main()
