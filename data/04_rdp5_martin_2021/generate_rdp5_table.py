"""
simulations/04_rdp5_martin_2021/generate_rdp5_table.py
======================================================
Generates LaTeX table tab_rdp5_benchmark.tex replicating Martin et al. (2021).
Adheres strictly to style guidelines: ZERO BOLD in prose and table.
"""

from pathlib import Path
import pandas as pd


def generate_latex_table(csv_path: Path, out_tex: Path):
    df = pd.read_csv(csv_path)

    tex = []
    tex.append("% LaTeX table generated from Martin et al. (2021) RDP5 benchmark suite")
    tex.append("% Strict compliance: zero bold formatting")
    tex.append(r"\begin{table*}[t!]")
    tex.append(r"\centering")
    tex.append(r"\small")
    tex.append(r"\caption{Computational scaling and runtime acceleration across the canonical RDP5 empirical benchmark suite. "
               r"Published runtimes for RDP4 and RDP5 are digitized from Martin et al.\ (2021) Supplementary Table S1 on an eight-core workstation. "
               r"\RhizAeon runtimes reflect wall-clock execution under the calibrated standard workflow on a standard laptop.}")
    tex.append(r"\label{tab:rdp5_suite_runtimes}")
    tex.append(r"{\setlength{\tabcolsep}{4.5pt}")
    tex.append(r"\begin{tabular*}{\textwidth}{@{\extracolsep{\fill}}l c c r r r r r@{}}")
    tex.append(r"\toprule")
    tex.append(r"Dataset & Taxa & Length (nt) & RDP4 & RDP5 & \RhizAeon & Speedup vs RDP5 & Speedup vs RDP4 \\")
    tex.append(r"\midrule")

    for _, row in df.iterrows():
        name = row["dataset_name"]
        taxa = int(row["taxa"])
        length = f"{int(row['length_nt']):,} nt"

        rdp4_str = f"{row['rdp4_published_sec']:.1f} s" if row['rdp4_published_sec'] < 60 else (
            f"{row['rdp4_published_sec']/60.0:.2f} m" if row['rdp4_published_sec'] < 3600 else f"{row['rdp4_published_sec']/3600.0:.2f} h"
        ) if pd.notna(row['rdp4_published_sec']) else "--"

        rdp5_str = f"{row['rdp5_published_sec']:.1f} s" if row['rdp5_published_sec'] < 60 else (
            f"{row['rdp5_published_sec']/60.0:.2f} m" if row['rdp5_published_sec'] < 3600 else f"{row['rdp5_published_sec']/3600.0:.2f} h"
        ) if pd.notna(row['rdp5_published_sec']) else "--"

        rhiz_t = row["rhiz_total_sec"]
        rhiz_str = f"{rhiz_t:.3f} s" if rhiz_t < 60 else f"{rhiz_t/60.0:.2f} m"

        sp5_str = f"{row['speedup_vs_rdp5']:.1f}$\\times$" if pd.notna(row['speedup_vs_rdp5']) else "Sub-second"
        sp4_str = f"{row['speedup_vs_rdp4']:.1f}$\\times$" if pd.notna(row['speedup_vs_rdp4']) else "Sub-second"

        tex.append(f"{name} & {taxa} & {length} & {rdp4_str} & {rdp5_str} & {rhiz_str} & {sp5_str} & {sp4_str} \\\\")

    tex.append(r"\bottomrule")
    tex.append(r"\end{tabular*}}")
    tex.append(r"\end{table*}")

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tex, "w") as f:
        f.write("\n".join(tex) + "\n")
    print(f"Saved LaTeX table to {out_tex}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    r_csv = base / "rdp5_benchmark_results.csv"
    tbl_path = Path(__file__).resolve().parent.parent.parent / "paper" / "tables" / "tab_rdp5_suite_runtimes.tex"
    if r_csv.exists():
        generate_latex_table(r_csv, tbl_path)
    else:
        print("CSV not found.")
