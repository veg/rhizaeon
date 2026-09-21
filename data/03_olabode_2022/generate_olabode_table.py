"""
simulations/03_olabode_2022/generate_olabode_table.py
=====================================================
Generates LaTeX table tab_olabode_benchmark.tex replicating Olabode et al. (2022).
Adheres strictly to style guidelines: ZERO BOLD in prose and table.
"""

from pathlib import Path
import numpy as np
import pandas as pd


def generate_latex_table(raw_csv: Path, summary_csv: Path, scaling_csv: Path, out_tex: Path):
    df_raw = pd.read_csv(raw_csv)
    df_sum = pd.read_csv(summary_csv)
    df_scale = pd.read_csv(scaling_csv) if scaling_csv.exists() else pd.DataFrame()

    tex = []
    tex.append("% LaTeX table generated from Olabode et al. (2022) replication benchmark")
    tex.append("% Strict compliance: zero bold formatting")
    tex.append(r"\begin{table}[htbp]")
    tex.append(r"\centering")
    tex.append(r"\small")
    tex.append(r"\caption{Performance of RhizAeon and 3SEQ across the Olabode et al. (2022) continuous-time HIV-1 benchmark suite. "
               r"Alignments span 9,000 nt under realistic continuous-time substitution histories with Gamma rate heterogeneity "
               r"($\alpha=1.5, \beta=3.0, \kappa=8.0$, tree length = 3.22 substitutions/codon). Partition misclassification error (\%) is evaluated "
               r"under optimal Hungarian bipartite matching against ground-truth phylogenetic partitions. Runtime denotes mean wall-clock latency per alignment.}")
    tex.append(r"\label{tab:olabode_benchmark}")
    tex.append(r"\begin{tabular}{llrrrrrrr}")
    tex.append(r"\toprule")
    tex.append(r"Experiment & Breakpoints & $N$ & Reps & Det. (\%) & Median Error (\%) & $Q_1$--$Q_3$ (\%) & Spatial (nt) & Runtime \\")
    tex.append(r"\midrule")

    # 1. Experiment 1: Post-hoc reference benchmark
    tex.append(r"\multicolumn{9}{l}{\textit{Experiment 1: Reference Genome Benchmark (COMET Design, $N=16$)}} \\")
    e1 = df_sum[df_sum["experiment"] == "exp1_posthoc"].sort_values("num_true_bps")
    for _, row in e1.iterrows():
        bp_label = f"{int(row['num_true_bps'])} Breakpoint" if row['num_true_bps'] == 1 else f"{int(row['num_true_bps'])} Breakpoints"
        iqr_str = f"[{row['rhiz_q25_error_pct']:.1f}, {row['rhiz_q75_error_pct']:.1f}]"
        # get mean spatial error from raw
        sub_raw = df_raw[(df_raw["experiment"] == "exp1_posthoc") & (df_raw["num_true_bps"] == row["num_true_bps"])]
        spat_err = sub_raw["rhiz_mean_spatial_error_nt"].dropna().mean()
        spat_str = f"{spat_err:.1f}" if not np.isnan(spat_err) else "--"
        time_str = f"{row['rhiz_mean_time_ms']:.1f} ms" if row['rhiz_mean_time_ms'] < 1000 else f"{row['rhiz_mean_time_ms']/1000.0:.2f} s"
        tex.append(f"Post-hoc & {bp_label} & {int(row['num_taxa'])} & {int(row['n_reps'])} & "
                   f"{row['rhiz_detection_rate']:.1f} & {row['rhiz_median_error_pct']:.2f} & {iqr_str} & {spat_str} & {time_str} \\\\")
    tex.append(r"\midrule")

    # 2. Experiment 2: Continuous time-scaled simulations with Gamma SRV
    tex.append(r"\multicolumn{9}{l}{\textit{Experiment 2: Continuous-Time HIV-1 Simulations with Gamma SRV ($N=37$)}} \\")
    e2 = df_sum[df_sum["experiment"] == "exp2_continuous"].sort_values("num_true_bps")
    for _, row in e2.iterrows():
        bp_label = f"{int(row['num_true_bps'])} Breakpoint" if row['num_true_bps'] == 1 else f"{int(row['num_true_bps'])} Breakpoints"
        iqr_str = f"[{row['rhiz_q25_error_pct']:.1f}, {row['rhiz_q75_error_pct']:.1f}]"
        sub_raw = df_raw[(df_raw["experiment"] == "exp2_continuous") & (df_raw["num_true_bps"] == row["num_true_bps"])]
        spat_err = sub_raw["rhiz_mean_spatial_error_nt"].dropna().mean()
        spat_str = f"{spat_err:.1f}" if not np.isnan(spat_err) else "--"
        time_str = f"{row['rhiz_mean_time_ms']:.1f} ms" if row['rhiz_mean_time_ms'] < 1000 else f"{row['rhiz_mean_time_ms']/1000.0:.2f} s"
        tex.append(f"Continuous & {bp_label} & {int(row['num_taxa'])} & {int(row['n_reps'])} & "
                   f"{row['rhiz_detection_rate']:.1f} & {row['rhiz_median_error_pct']:.2f} & {iqr_str} & {spat_str} & {time_str} \\\\")
    tex.append(r"\midrule")

    # 3. Experiment 3: Scaling Cohorts
    if not df_scale.empty:
        tex.append(r"\multicolumn{9}{l}{\textit{Experiment 3: Computational Scaling Across Taxa ($L=9{,}000$ nt)}} \\")
        for _, row in df_scale.iterrows():
            sub_raw = df_raw[(df_raw["experiment"] == "exp3_scaling") & (df_raw["num_taxa"] == row["num_taxa"])]
            det_rate = sub_raw["rhiz_detected"].mean() * 100.0
            med_err = sub_raw["rhiz_error_pct"].median()
            q25 = sub_raw["rhiz_error_pct"].quantile(0.25)
            q75 = sub_raw["rhiz_error_pct"].quantile(0.75)
            iqr_str = f"[{q25:.1f}, {q75:.1f}]"
            spat_str = "--"
            time_str = f"{row['rhiz_time_sec']:.2f} s"
            tex.append(f"Scaling & 1--3 BP & {int(row['num_taxa'])} & {int(row['n_reps'])} & "
                       f"{det_rate:.1f} & {med_err:.2f} & {iqr_str} & {spat_str} & {time_str} \\\\")

    tex.append(r"\bottomrule")
    tex.append(r"\end{tabular}")
    tex.append(r"\end{table}")

    out_tex.parent.mkdir(parents=True, exist_ok=True)
    with open(out_tex, "w") as f:
        f.write("\n".join(tex) + "\n")
    print(f"Saved LaTeX table to {out_tex}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent
    r_csv = base / "olabode_raw_results.csv"
    s_csv = base / "olabode_summary.csv"
    sc_csv = base / "olabode_scaling.csv"
    tbl_path = Path(__file__).resolve().parent.parent.parent / "paper" / "tables" / "tab_olabode_benchmark.tex"
    if r_csv.exists() and s_csv.exists():
        generate_latex_table(r_csv, s_csv, sc_csv, tbl_path)
    else:
        print("Summary CSV not yet generated.")
