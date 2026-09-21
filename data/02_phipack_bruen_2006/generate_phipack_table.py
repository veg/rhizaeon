"""
simulations/02_phipack_bruen_2006/generate_phipack_table.py
============================================================
Generates LaTeX table tab_phipack_benchmark.tex replicating Bruen et al. (2006).
Adheres strictly to style guidelines: ZERO BOLD in prose and table.
"""

from pathlib import Path
import pandas as pd


def generate_latex_table(summary_csv: Path, out_tex: Path):
    df = pd.read_csv(summary_csv)

    tex = []
    tex.append("% LaTeX table generated from Bruen et al. (2006) replication benchmark")
    tex.append("% Strict compliance: zero bold formatting")
    tex.append(r"\begin{table}[htbp]")
    tex.append(r"\centering")
    tex.append(r"\small")
    tex.append(r"\caption{Performance of RhizAeon, PHI, Max$\chi^2$, NSS, and 3SEQ across the Bruen et al. (2006) coalescent benchmark suite. "
               r"Alignments span 1,000 nt under neutral constant size, exponential population growth ($\beta$), Gamma rate heterogeneity ($\alpha$), "
               r"and varying recombination intensity ($\rho = 4N_0 r$) across sample sizes $n \in \{10, 50\}$. False positive rates (FPR) and detection power "
               r"are evaluated across 50 replicates per parameter combination.}")
    tex.append(r"\label{tab:phipack_benchmark}")
    tex.append(r"\begin{tabular}{lllrrrrrr}")
    tex.append(r"\toprule")
    tex.append(r"Regime & $n$ & Parameters & Inf. Sites & RhizAeon (\%) & PHI (\%) & Max$\chi^2$ (\%) & NSS (\%) & 3SEQ (\%) \\")
    tex.append(r"\midrule")

    # 1. Null Constant
    tex.append(r"\multicolumn{9}{l}{\textit{Neutral Constant Population Size ($\rho = 0$)}} \\")
    c_df = df[df["category"] == "null_constant"]
    for _, row in c_df.iterrows():
        p_str = f"$\\theta={int(row['theta'])}$"
        tex.append(f"Neutral & {int(row['n_taxa'])} & {p_str} & {row['mean_informative_sites']:.1f} & "
                   f"{row['rhiz_rate']:.1f} & {row['phi_rate']:.1f} & {row['maxchi_rate']:.1f} & {row['nss_rate']:.1f} & {row['three_rate']:.1f} \\\\")
    tex.append(r"\midrule")

    # 2. Null Exponential Growth
    tex.append(r"\multicolumn{9}{l}{\textit{Exponential Population Expansion ($\rho = 0$, Star-like Genealogies)}} \\")
    g_df = df[df["category"] == "null_growth"].sort_values(["n_taxa", "beta", "theta"])
    for _, row in g_df.iterrows():
        p_str = f"$\\beta={int(row['beta'])}, \\theta={int(row['theta'])}$"
        tex.append(f"Exp. Growth & {int(row['n_taxa'])} & {p_str} & {row['mean_informative_sites']:.1f} & "
                   f"{row['rhiz_rate']:.1f} & {row['phi_rate']:.1f} & {row['maxchi_rate']:.1f} & {row['nss_rate']:.1f} & {row['three_rate']:.1f} \\\\")
    tex.append(r"\midrule")

    # 3. Null Gamma Rate Heterogeneity
    tex.append(r"\multicolumn{9}{l}{\textit{Gamma Site Rate Variation ($\rho = 0$)}} \\")
    gam_df = df[df["category"] == "null_gamma"].sort_values(["n_taxa", "alpha", "theta"])
    for _, row in gam_df.iterrows():
        p_str = f"$\\alpha={row['alpha']}, \\theta={int(row['theta'])}$"
        tex.append(f"Gamma & {int(row['n_taxa'])} & {p_str} & {row['mean_informative_sites']:.1f} & "
                   f"{row['rhiz_rate']:.1f} & {row['phi_rate']:.1f} & {row['maxchi_rate']:.1f} & {row['nss_rate']:.1f} & {row['three_rate']:.1f} \\\\")
    tex.append(r"\midrule")

    # 4. Recombination Power
    tex.append(r"\multicolumn{9}{l}{\textit{Recombination Power ($\theta = 20, \rho > 0$)}} \\")
    pow_df = df[(df["category"] == "power") & (df["theta"] == 20.0) & (df["rho"].isin([4.0, 16.0, 64.0]))].sort_values(["n_taxa", "rho"])
    for _, row in pow_df.iterrows():
        p_str = f"$\\rho={int(row['rho'])}$"
        tex.append(f"Power & {int(row['n_taxa'])} & {p_str} & {row['mean_informative_sites']:.1f} & "
                   f"{row['rhiz_rate']:.1f} & {row['phi_rate']:.1f} & {row['maxchi_rate']:.1f} & {row['nss_rate']:.1f} & {row['three_rate']:.1f} \\\\")

    tex.append(r"\bottomrule")
    tex.append(r"\end{tabular}")
    tex.append(r"\end{table}")

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tex, "w") as f:
        f.write("\n".join(tex) + "\n")
    print(f"Saved LaTeX table to {out_tex}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    sum_csv = base / "phipack_summary.csv"
    tbl_path = Path(__file__).resolve().parent.parent.parent / "paper" / "tables" / "tab_phipack_benchmark.tex"
    if sum_csv.exists():
        generate_latex_table(sum_csv, tbl_path)
    else:
        print("Summary CSV not yet generated.")
