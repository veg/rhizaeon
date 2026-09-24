#!/usr/bin/env python3
"""
build_intuition_page.py
========================
Compiles docs/rhizaeon_overview_and_intuition.md into a responsive,
publication-quality HTML dossier (intuition.html) in rhizaeon_bench/
adhering to the BRC-Analytics / RhizAeon design system.
"""

import os
import re
import subprocess
from pathlib import Path

BENCH_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination/rhizaeon_bench")
ROOT_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination")
MD_SOURCE = ROOT_DIR / "docs" / "rhizaeon_overview_and_intuition.md"
HTML_OUTPUT = BENCH_DIR / "intuition.html"

def main():
    print(f"[*] Reading source Markdown from {MD_SOURCE}...")
    with open(MD_SOURCE, "r", encoding="utf-8") as f:
        md_text = f.read()

    # Update relative image paths to point to assets/figures/
    md_text_web = re.sub(r'\.\./paper/figures/([^)]+)', r'assets/figures/\1', md_text)

    temp_md = BENCH_DIR / "scratch_intuition.md"
    with open(temp_md, "w", encoding="utf-8") as f:
        f.write(md_text_web)

    print("[*] Running Pandoc to generate HTML body...")
    cmd = [
        "pandoc",
        str(temp_md),
        "-f", "markdown+pipe_tables+backtick_code_blocks+tex_math_dollars+raw_html",
        "-t", "html",
        "--mathjax"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("[!] Pandoc error:", res.stderr)
        return

    html_body = res.stdout

    # Assemble complete page with BRC-Analytics design system
    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>RhizAeon Methodological &amp; Geometric Intuition Guide | BRC-Analytics</title>
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=Roboto+Mono:wght@400;500;600&display=swap" rel="stylesheet">
  <link rel="stylesheet" href="assets/css/style.css">
  <script>
    window.MathJax = {{
      tex: {{
        inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
        displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']],
        processEscapes: true
      }},
      options: {{
        skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
      }}
    }};
  </script>
  <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
  <style>
    .intuition-hero {{
      background: linear-gradient(135deg, #1e3a8a 0%, #0369a1 100%);
      color: #ffffff;
      padding: 3rem 0 2.5rem 0;
      border-bottom: 1px solid #0284c7;
    }}
    .intuition-eyebrow {{
      font-size: 0.8rem;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      color: #bae6fd;
      margin-bottom: 0.5rem;
    }}
    .intuition-title {{
      font-size: 2.25rem;
      font-weight: 800;
      line-height: 1.2;
      margin-bottom: 0.8rem;
      color: #ffffff;
      letter-spacing: -0.02em;
    }}
    .intuition-sub {{
      font-size: 1.05rem;
      color: #e0f2fe;
      line-height: 1.55;
      max-width: 880px;
      margin-bottom: 1.5rem;
    }}
    .intuition-action-bar {{
      display: flex;
      flex-wrap: wrap;
      gap: 1rem;
      align-items: center;
    }}
    .btn-pdf-download {{
      background: #ffffff;
      color: #0369a1;
      font-weight: 700;
      font-size: 0.9rem;
      padding: 0.65rem 1.25rem;
      border-radius: 6px;
      display: inline-flex;
      align-items: center;
      gap: 0.5rem;
      box-shadow: 0 2px 6px rgba(0,0,0,0.15);
      transition: all 0.15s ease;
      text-decoration: none;
    }}
    .btn-pdf-download:hover {{
      background: #f0f9ff;
      transform: translateY(-1px);
      text-decoration: none;
    }}
    .btn-outline {{
      background: rgba(255,255,255,0.12);
      border: 1px solid rgba(255,255,255,0.4);
      color: #ffffff;
      font-weight: 600;
      font-size: 0.9rem;
      padding: 0.65rem 1.15rem;
      border-radius: 6px;
      display: inline-flex;
      align-items: center;
      gap: 0.4rem;
      transition: all 0.15s ease;
      text-decoration: none;
    }}
    .btn-outline:hover {{
      background: rgba(255,255,255,0.22);
      color: #ffffff;
      text-decoration: none;
    }}
    .doc-layout {{
      display: grid;
      grid-template-columns: 280px 1fr;
      gap: 2.5rem;
      padding: 2.5rem 0;
    }}
    @media (max-width: 992px) {{
      .doc-layout {{
        grid-template-columns: 1fr;
      }}
    }}
    .doc-sidebar {{
      position: sticky;
      top: 5rem;
      align-self: start;
      background: #ffffff;
      border: 1px solid #e2e8f0;
      border-radius: 8px;
      padding: 1.25rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
      font-size: 0.85rem;
    }}
    .doc-sidebar-title {{
      font-weight: 700;
      color: #0f172a;
      text-transform: uppercase;
      font-size: 0.72rem;
      letter-spacing: 0.06em;
      margin-bottom: 0.75rem;
      border-bottom: 1px solid #e2e8f0;
      padding-bottom: 0.4rem;
    }}
    .doc-sidebar ul {{
      list-style: none;
      padding: 0;
      margin: 0;
    }}
    .doc-sidebar li {{
      margin-bottom: 0.5rem;
    }}
    .doc-sidebar a {{
      color: #475569;
      text-decoration: none;
      display: block;
      padding: 0.2rem 0;
      line-height: 1.35;
      transition: color 0.15s ease;
    }}
    .doc-sidebar a:hover {{
      color: #0284c7;
      text-decoration: underline;
    }}
    .doc-content {{
      background: #ffffff;
      border: 1px solid #e2e8f0;
      border-radius: 8px;
      padding: 2.5rem 3rem;
      box-shadow: 0 1px 3px rgba(0,0,0,0.05);
      line-height: 1.68;
      font-size: 10pt;
      color: #1e293b;
    }}
    .doc-content h1 {{
      font-size: 2rem;
      color: #0f172a;
      margin-top: 0;
      margin-bottom: 1.2rem;
      padding-bottom: 0.6rem;
      border-bottom: 2px solid #0284c7;
      font-weight: 800;
    }}
    .doc-content h2 {{
      font-size: 1.4rem;
      color: #1e3a8a;
      margin-top: 2.2rem;
      margin-bottom: 0.8rem;
      padding-bottom: 0.4rem;
      border-bottom: 1px solid #cbd5e1;
      font-weight: 700;
    }}
    .doc-content h3 {{
      font-size: 1.15rem;
      color: #0369a1;
      margin-top: 1.6rem;
      margin-bottom: 0.6rem;
      font-weight: 600;
    }}
    .doc-content h4 {{
      font-size: 1.02rem;
      color: #334155;
      margin-top: 1.2rem;
      margin-bottom: 0.4rem;
      font-weight: 600;
    }}
    .doc-content p {{
      margin-bottom: 1.1rem;
      text-align: justify;
    }}
    .doc-content ul, .doc-content ol {{
      margin: 1rem 0 1.2rem 1.8rem;
    }}
    .doc-content li {{
      margin-bottom: 0.45rem;
    }}
    .doc-content pre {{
      background: #f8fafc;
      border: 1px solid #cbd5e1;
      border-left: 4px solid #0284c7;
      border-radius: 4px;
      padding: 1rem 1.25rem;
      font-family: var(--font-mono);
      font-size: 8.2pt;
      line-height: 1.35;
      overflow-x: auto;
      margin: 1.4rem 0;
    }}
    .doc-content code {{
      font-family: var(--font-mono);
      background: #f1f5f9;
      color: #0f172a;
      padding: 2px 6px;
      border-radius: 3px;
      font-size: 8.8pt;
      border: 1px solid #e2e8f0;
    }}
    .doc-content pre code {{
      background: transparent;
      padding: 0;
      border: none;
      font-size: inherit;
    }}
    .doc-content table {{
      width: 100%;
      border-collapse: collapse;
      margin: 1.5rem 0;
      font-size: 8.8pt;
    }}
    .doc-content th {{
      background: #f1f5f9;
      color: #0f172a;
      font-weight: 700;
      text-align: left;
      padding: 8px 12px;
      border-top: 1.5px solid #64748b;
      border-bottom: 1.5px solid #64748b;
    }}
    .doc-content td {{
      padding: 8px 12px;
      border-bottom: 1px solid #e2e8f0;
      vertical-align: top;
    }}
    .doc-content tr:nth-child(even) td {{
      background: #f8fafc;
    }}
    .doc-content figure {{
      margin: 2rem 0;
      text-align: center;
      background: #f8fafc;
      border: 1px solid #e2e8f0;
      border-radius: 6px;
      padding: 1.25rem;
    }}
    .doc-content img {{
      max-width: 100%;
      height: auto;
      border-radius: 4px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      margin-bottom: 0.8rem;
    }}
    .doc-content figcaption {{
      font-size: 8.5pt;
      color: #475569;
      line-height: 1.45;
      text-align: justify;
      margin-top: 0.5rem;
    }}
    .doc-content blockquote {{
      border-left: 4px solid #0284c7;
      background: #f0f9ff;
      margin: 1.4rem 0;
      padding: 1rem 1.25rem;
      border-radius: 0 6px 6px 0;
      color: #0369a1;
    }}
    .doc-content blockquote p {{
      margin: 0;
    }}
  </style>
</head>
<body>

  <!-- Site Navigation -->
  <header class="rams-header">
    <div class="rams-container rams-nav-flex">
      <a href="index.html" class="rams-brand">
        <span class="rams-brand-name">RhizAeon</span>
        <span class="rams-brand-sep">/</span>
        <span class="rams-brand-desc">Method Intuition &amp; Architecture Guide</span>
      </a>
      <nav class="rams-nav">
        <a href="index.html">&larr; Benchmark Compendium</a>
        <a href="rhizaeon_overview_and_intuition.pdf" target="_blank" rel="noopener" class="rams-nav-chip">Download PDF (26 pp) &nearr;</a>
        <a href="https://brc-analytics.org" target="_blank" rel="noopener">BRC-Analytics &nearr;</a>
        <a href="AGENT.MD">AGENT.md</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">GitHub &nearr;</a>
      </nav>
    </div>
  </header>

  <!-- Hero Banner -->
  <section class="intuition-hero">
    <div class="rams-container">
      <div class="intuition-eyebrow">THEORETICAL TREATISE &bull; METHODOLOGICAL INTUITION</div>
      <h1 class="intuition-title">How RhizAeon Works: Coordinate Manifolds, Dynamic Trajectories &amp; Dual-Tier Architecture</h1>
      <p class="intuition-sub">
        A comprehensive illustrated guide to detecting viral and bacterial recombination directly from sequence coordinate manifolds—without guide trees, phylogenetic heuristics, or Markov chain Monte Carlo simulations.
      </p>
      <div class="intuition-action-bar">
        <a href="rhizaeon_overview_and_intuition.pdf" target="_blank" rel="noopener" class="btn-pdf-download">
          <svg width="18" height="18" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path>
            <polyline points="7 10 12 15 17 10"></polyline>
            <line x1="12" y1="15" x2="12" y2="3"></line>
          </svg>
          <span>Download Publication PDF (26 Pages, 3.5 MB)</span>
        </a>
        <a href="index.html#benchmarks" class="btn-outline">
          <span>Explore 11 Benchmark Suites &rarr;</span>
        </a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener" class="btn-outline">
          <span>Source Code on GitHub &nearr;</span>
        </a>
      </div>
    </div>
  </section>

  <!-- Main Content Layout -->
  <div class="rams-container doc-layout">
    
    <!-- Sidebar Navigation -->
    <aside class="doc-sidebar">
      <div class="doc-sidebar-title">Guide Sections</div>
      <ul>
        <li><a href="#executive-intuition">1. Executive Intuition</a></li>
        <li><a href="#why-phylogenetics-fails-at-planetary-scale">2. The Phylogenetic Scaling Crisis</a></li>
        <li><a href="#tier-1-geometric-manifold-engine">3. Tier 1 Geometric Engine</a></li>
        <li><a href="#the-three-core-mathematical-operations">4. Three Core Mathematical Operations</a></li>
        <li><a href="#the-canonical-walkthrough-hiv-1-kal153">5. Walkthrough: HIV-1 KAL153</a></li>
        <li><a href="#tier-1-breakdown-triggers-failure-mode-analysis">6. Tier 1 Breakdown Triggers</a></li>
        <li><a href="#tier-2-attention-and-conformal-deep-resolution">7. Tier 2 Attention &amp; Conformal Engine</a></li>
        <li><a href="#empirical-validation-across-sars-cov-2-spike">8. SARS-CoV-2 Spike Validation</a></li>
        <li><a href="#computational-complexity-scaling-laws">9. Scaling Laws &amp; Complexity</a></li>
        <li><a href="#theoretical-synthesis-the-new-recombination-paradigm">10. Theoretical Synthesis</a></li>
      </ul>
      <div style="margin-top: 1.5rem; pt: 1rem; border-top: 1px solid #e2e8f0;">
        <a href="rhizaeon_overview_and_intuition.pdf" target="_blank" rel="noopener" style="font-weight: 700; color: #0284c7; display: flex; align-items: center; gap: 0.35rem;">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round">
            <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"></path>
            <polyline points="14 2 14 8 20 8"></polyline>
          </svg>
          <span>Full 26-Page PDF &nearr;</span>
        </a>
      </div>
    </aside>

    <!-- Document Content -->
    <article class="doc-content">
      {html_body}
    </article>

  </div>

  <!-- Site Footer -->
  <footer class="rams-footer" style="background: #0f172a; color: #94a3b8; padding: 3rem 0; margin-top: 4rem; border-top: 1px solid #1e293b;">
    <div class="rams-container" style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 1.5rem;">
      <div>
        <div style="color: #ffffff; font-weight: 700; font-size: 1.1rem; margin-bottom: 0.3rem;">RhizAeon Continuous Coordinate Recombination Framework</div>
        <p style="font-size: 0.85rem; margin: 0; color: #64748b;">
          Supported by the National Institute of Allergy and Infectious Diseases (NIAID) under BRC-Analytics.
        </p>
      </div>
      <div style="display: flex; gap: 1.5rem; font-size: 0.85rem;">
        <a href="index.html" style="color: #cbd5e1;">Home</a>
        <a href="rhizaeon_overview_and_intuition.pdf" target="_blank" rel="noopener" style="color: #cbd5e1;">PDF Guide</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener" style="color: #cbd5e1;">GitHub Repository</a>
      </div>
    </div>
  </footer>

</body>
</html>
"""

    with open(HTML_OUTPUT, "w", encoding="utf-8") as f:
        f.write(full_html)

    if temp_md.exists():
        temp_md.unlink()

    print(f"[✓] Successfully compiled {HTML_OUTPUT}")
    print(f"[✓] File size: {os.path.getsize(HTML_OUTPUT):,} bytes")

if __name__ == "__main__":
    main()
