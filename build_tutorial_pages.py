#!/usr/bin/env python3
"""
build_tutorial_pages.py
========================
Compiles the four tutorial tracks from docs/ into standalone responsive HTML
pages adhering to the BRC-Analytics / RhizAeon design system:
  - track1_visual_quickstart.html
  - track2_practitioner_guide.html
  - track3_algorithmic_derivation.html
  - track4_conceptual_monograph.html
  - tutorials.html (Master Diataxis index)
"""

import os
import re
import subprocess
from pathlib import Path

BENCH_DIR = Path("/Users/sergei/Projects/TOGA_MEME/recombination/rhizaeon_bench")
DOCS_DIR = BENCH_DIR / "docs"

TRACKS = [
    {
        "id": "track1",
        "src": DOCS_DIR / "track1_visual_quickstart.md",
        "dst": BENCH_DIR / "track1_visual_quickstart.html",
        "title": "Track 1: Fast-Track Worked-Example (Show-First)",
        "eyebrow": "TUTORIAL TRACK 1 • SHOW-FIRST QUICKSTART",
        "desc": "Witness RhizAeon in flight on canonical HIV-1 KAL153 in under 2 minutes. Zero abstract derivations before the primary visual reveal."
    },
    {
        "id": "track2",
        "src": DOCS_DIR / "track2_practitioner_guide.md",
        "dst": BENCH_DIR / "track2_practitioner_guide.html",
        "title": "Track 2: Diagnostic & Forensic Workflow (Practitioner How-To)",
        "eyebrow": "TUTORIAL TRACK 2 • PRACTITIONER HOW-TO",
        "desc": "A clinical diagnostic guide for real-world genomic data: QC gates, parameter trade-off matrices, troubleshooting triage, and alignment disassembly."
    },
    {
        "id": "track3",
        "src": DOCS_DIR / "track3_algorithmic_derivation.md",
        "dst": BENCH_DIR / "track3_algorithmic_derivation.html",
        "title": "Track 3: Deductive / Algorithmic Derivation (Mathematical Foundations)",
        "eyebrow": "TUTORIAL TRACK 3 • MATHEMATICAL DERIVATION",
        "desc": "Rigorous axiomatic problem formulation, Gram matrix eigendecomposition, fixed projection operator proofs, and O(N² L) complexity bounds."
    },
    {
        "id": "track4",
        "src": DOCS_DIR / "track4_conceptual_monograph.md",
        "dst": BENCH_DIR / "track4_conceptual_monograph.html",
        "title": "Track 4: Conceptual Monograph & Operational Boundaries",
        "eyebrow": "TUTORIAL TRACK 4 • ARCHITECTURAL DEEP-DIVE",
        "desc": "The 1,050 simulation envelope, Frame Rigidity collapse, Graph Laplacian hand-off, circulating recombinant forms (CRFs), and the Evolutionary Sieve."
    },
    {
        "id": "index",
        "src": DOCS_DIR / "README.md",
        "dst": BENCH_DIR / "tutorials.html",
        "title": "RhizAeon Multi-Tier Tutorial Suite | BRC-Analytics",
        "eyebrow": "DOCUMENTATION SUITE • DIÁTAXIS ARCHITECTURE",
        "desc": "A synchronized suite of four tutorial tracks tailored to different cognitive modes, user proficiencies, and learning objectives."
    }
]


def render_html_page(track_info):
    src_file = track_info["src"]
    dst_file = track_info["dst"]

    print(f"[*] Compiling {src_file.name} -> {dst_file.name}...")
    with open(src_file, "r", encoding="utf-8") as f:
        md_text = f.read()

    # Rewrite image links to point to assets/figures/
    md_text_web = re.sub(r'\.\./paper/figures/([^)]+)', r'assets/figures/\1', md_text)
    md_text_web = re.sub(r'track1_visual_quickstart\.md', r'track1_visual_quickstart.html', md_text_web)
    md_text_web = re.sub(r'track2_practitioner_guide\.md', r'track2_practitioner_guide.html', md_text_web)
    md_text_web = re.sub(r'track3_algorithmic_derivation\.md', r'track3_algorithmic_derivation.html', md_text_web)
    md_text_web = re.sub(r'track4_conceptual_monograph\.md', r'track4_conceptual_monograph.html', md_text_web)

    temp_md = BENCH_DIR / f"temp_{track_info['id']}.md"
    with open(temp_md, "w", encoding="utf-8") as f:
        f.write(md_text_web)

    cmd = [
        "pandoc",
        str(temp_md),
        "-f", "markdown+pipe_tables+backtick_code_blocks+tex_math_dollars+raw_html",
        "-t", "html",
        "--mathjax"
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"[!] Pandoc error on {src_file.name}:", res.stderr)
        return

    html_body = res.stdout

    full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{track_info['title']}</title>
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
    .track-hero {{
      background: linear-gradient(135deg, #0f172a 0%, #0369a1 100%);
      color: #ffffff;
      padding: 3rem 0 2.5rem 0;
      border-bottom: 1px solid #0284c7;
    }}
    .track-eyebrow {{
      font-size: 0.8rem;
      font-weight: 700;
      letter-spacing: 0.08em;
      text-transform: uppercase;
      color: #7dd3fc;
      margin-bottom: 0.5rem;
    }}
    .track-title {{
      font-size: 2.2rem;
      font-weight: 800;
      line-height: 1.2;
      margin-bottom: 0.8rem;
      color: #ffffff;
      letter-spacing: -0.02em;
    }}
    .track-sub {{
      font-size: 1.05rem;
      color: #e0f2fe;
      line-height: 1.55;
      max-width: 900px;
      margin-bottom: 1.5rem;
    }}
    .track-nav-tabs {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      margin-top: 1rem;
    }}
    .track-tab {{
      padding: 0.45rem 0.9rem;
      border-radius: 6px;
      font-size: 0.82rem;
      font-weight: 600;
      text-decoration: none;
      transition: all 0.15s ease;
      background: rgba(255,255,255,0.12);
      color: #e0f2fe;
      border: 1px solid rgba(255,255,255,0.25);
    }}
    .track-tab:hover {{
      background: rgba(255,255,255,0.25);
      color: #ffffff;
      text-decoration: none;
    }}
    .track-tab.active {{
      background: #38bdf8;
      color: #0f172a;
      border-color: #38bdf8;
      font-weight: 700;
    }}
    .doc-layout {{
      display: grid;
      grid-template-columns: 260px 1fr;
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
      line-height: 1.7;
      color: #1e293b;
      min-width: 0;
    }}
    .doc-content h1, .doc-content h2, .doc-content h3, .doc-content h4 {{
      color: #0f172a;
      font-weight: 700;
      margin-top: 2rem;
      margin-bottom: 0.75rem;
      letter-spacing: -0.015em;
    }}
    .doc-content h1 {{ font-size: 1.85rem; border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; }}
    .doc-content h2 {{ font-size: 1.45rem; border-bottom: 1px solid #e2e8f0; padding-bottom: 0.35rem; }}
    .doc-content h3 {{ font-size: 1.15rem; }}
    .doc-content pre {{
      background: #0f172a;
      color: #f8fafc;
      padding: 1.1rem 1.3rem;
      border-radius: 6px;
      overflow-x: auto;
      font-family: 'Roboto Mono', monospace;
      font-size: 0.85rem;
      line-height: 1.5;
      margin: 1.25rem 0;
    }}
    .doc-content code {{
      font-family: 'Roboto Mono', monospace;
      font-size: 0.88em;
      background: #f1f5f9;
      color: #0f172a;
      padding: 0.15em 0.35em;
      border-radius: 4px;
      border: 1px solid #e2e8f0;
    }}
    .doc-content pre code {{
      background: transparent;
      color: inherit;
      padding: 0;
      border: none;
    }}
    .doc-content table {{
      width: 100%;
      border-collapse: collapse;
      margin: 1.5rem 0;
      font-size: 0.88rem;
    }}
    .doc-content th, .doc-content td {{
      border: 1px solid #e2e8f0;
      padding: 0.65rem 0.9rem;
      text-align: left;
    }}
    .doc-content th {{
      background: #f8fafc;
      font-weight: 600;
      color: #0f172a;
    }}
    .doc-content img {{
      max-width: 100%;
      height: auto;
      border-radius: 6px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
      margin: 1.5rem 0;
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
        <span class="rams-brand-desc">Multi-Tier Tutorial Suite</span>
      </a>
      <nav class="rams-nav">
        <a href="index.html">&larr; Benchmark Compendium</a>
        <a href="tutorials.html" class="rams-nav-chip" style="background: #38bdf8; color: #0f172a; font-weight: 700;">Tutorial Suite Index</a>
        <a href="intuition.html">Monograph (26 pp)</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener">GitHub &nearr;</a>
      </nav>
    </div>
  </header>

  <!-- Hero Banner -->
  <section class="track-hero">
    <div class="rams-container">
      <div class="track-eyebrow">{track_info['eyebrow']}</div>
      <h1 class="track-title">{track_info['title']}</h1>
      <p class="track-sub">{track_info['desc']}</p>
      
      <!-- Track Switcher -->
      <div class="track-nav-tabs">
        <a href="tutorials.html" class="track-tab {'active' if track_info['id'] == 'index' else ''}">Overview Index</a>
        <a href="track1_visual_quickstart.html" class="track-tab {'active' if track_info['id'] == 'track1' else ''}">Track 1: Quickstart</a>
        <a href="track2_practitioner_guide.html" class="track-tab {'active' if track_info['id'] == 'track2' else ''}">Track 2: Practitioner</a>
        <a href="track3_algorithmic_derivation.html" class="track-tab {'active' if track_info['id'] == 'track3' else ''}">Track 3: Mathematical</a>
        <a href="track4_conceptual_monograph.html" class="track-tab {'active' if track_info['id'] == 'track4' else ''}">Track 4: Monograph</a>
      </div>
    </div>
  </section>

  <!-- Main Content Layout -->
  <div class="rams-container doc-layout">
    
    <!-- Sidebar Navigation -->
    <aside class="doc-sidebar">
      <div class="doc-sidebar-title">Tutorial Navigation</div>
      <ul>
        <li><a href="tutorials.html" style="font-weight: {'700' if track_info['id'] == 'index' else '400'}; color: {'#0284c7' if track_info['id'] == 'index' else '#475569'};">&bull; Suite Overview</a></li>
        <li><a href="track1_visual_quickstart.html" style="font-weight: {'700' if track_info['id'] == 'track1' else '400'}; color: {'#0284c7' if track_info['id'] == 'track1' else '#475569'};">&bull; Track 1: Quickstart</a></li>
        <li><a href="track2_practitioner_guide.html" style="font-weight: {'700' if track_info['id'] == 'track2' else '400'}; color: {'#0284c7' if track_info['id'] == 'track2' else '#475569'};">&bull; Track 2: Practitioner</a></li>
        <li><a href="track3_algorithmic_derivation.html" style="font-weight: {'700' if track_info['id'] == 'track3' else '400'}; color: {'#0284c7' if track_info['id'] == 'track3' else '#475569'};">&bull; Track 3: Mathematical</a></li>
        <li><a href="track4_conceptual_monograph.html" style="font-weight: {'700' if track_info['id'] == 'track4' else '400'}; color: {'#0284c7' if track_info['id'] == 'track4' else '#475569'};">&bull; Track 4: Monograph</a></li>
      </ul>
      <div style="margin-top: 1.5rem; padding-top: 1rem; border-top: 1px solid #e2e8f0;">
        <a href="intuition.html" style="font-weight: 600; color: #0284c7; display: block; margin-bottom: 0.4rem;">&rarr; 26-Page Monograph</a>
        <a href="rhizaeon_overview_and_intuition.pdf" target="_blank" rel="noopener" style="font-weight: 600; color: #0284c7; display: block;">&rarr; Publication PDF &nearr;</a>
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
        <a href="tutorials.html" style="color: #cbd5e1;">Tutorial Suite</a>
        <a href="intuition.html" style="color: #cbd5e1;">Monograph</a>
        <a href="https://github.com/veg/rhizaeon" target="_blank" rel="noopener" style="color: #cbd5e1;">GitHub Repository</a>
      </div>
    </div>
  </footer>

</body>
</html>
"""

    with open(dst_file, "w", encoding="utf-8") as f:
        f.write(full_html)

    if temp_md.exists():
        temp_md.unlink()

    print(f"[✓] Created {dst_file.name} ({os.path.getsize(dst_file):,} bytes)")


def main():
    print("[*] Generating all 5 tutorial web pages from Markdown...")
    for t in TRACKS:
        render_html_page(t)
    print("\n[✓] All tutorial pages compiled successfully!")


if __name__ == "__main__":
    main()
