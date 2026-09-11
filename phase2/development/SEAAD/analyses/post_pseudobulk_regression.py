#!/usr/bin/env python3
"""
Post-process disease pseudobulk regression results.
Aggregates cell-type results, computes FDR corrections, and generates volcano plots.
"""

import sys
import logging
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas import DataFrame, concat, read_csv
from statsmodels.stats.multitest import multipletests

# Configure logging
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run post-processing for disease pseudobulk regression analysis."
    )
    parser.add_argument(
        "--project",
        type=str,
        default="seaad_ec_multiome",
        help="Project name used for file prefixes.",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad",
        help="Base working directory.",
    )
    parser.add_argument(
        "--modality",
        type=str,
        default="rna",
        choices=["rna", "atac"],
        help="Data modality (rna or atac).",
    )
    parser.add_argument(
        "--target-variable",
        type=str,
        default="dx",
        help="The primary disease target variable column in covariates (default: 'dx').",
    )
    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
        choices=["ols", "glm", "glm_tweedie", "rlm", "wls", "vwrlm"],
        help="Regression method to use.",
    )
    parser.add_argument(
        "--min-rlm-effect",
        type=float,
        default=0.001,
        help="Minimum effect size for RLM significant results.",
    )
    parser.add_argument(
        "--no-volcano-per-celltype",
        action="store_false",
        dest="volcano_per_celltype",
        help="Disable generation of volcano plots per cell-type.",
    )
    parser.set_defaults(volcano_per_celltype=True)
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def compute_bh_fdr(
    df: DataFrame,
    alpha: float = 0.05,
    p_col: str = "p-value",
    method: str = "fdr_bh",
    verbose: bool = True,
) -> DataFrame:
    ret_df = df.copy()
    p_vals = ret_df[p_col].fillna(1)
    test_adjust = multipletests(np.array(p_vals), alpha=alpha, method=method)
    ret_df[method] = test_adjust[1]
    if verbose:
        logger.info(
            "Total significant after correction: %d",
            ret_df.loc[ret_df[method] < alpha].shape[0],
        )
    return ret_df


def volcano_plot(
    df: DataFrame,
    project: str,
    modality: str,
    regression_type: str,
    target_variable: str,
    figures_dir: Path,
    x_term: str = "log2fc",
    y_term: str = "p-value",
    alpha: float = 0.05,
    adj_p_col: str = "fdr_bh",
    title: str = None,
    filter_nseeff: bool = True,
    extreme_size: float = 10.0,
):
    plot_df = df.copy()
    plot_df = plot_df.reset_index(drop=True)

    if filter_nseeff:
        plot_df = plot_df.loc[
            (
                (-extreme_size < plot_df[x_term])
                & (plot_df[x_term] < extreme_size)
                & (~plot_df["z"].isna())
                | (plot_df[adj_p_col] < alpha)
            )
        ]

    plt.figure(figsize=(9, 9))

    plot_df[y_term] = plot_df[y_term].replace(0, np.finfo(float).eps)
    log_pvalue = -np.log10(plot_df[y_term])

    is_sig = plot_df[adj_p_col] < alpha

    sns.set_style("whitegrid")
    sns.scatterplot(
        x=plot_df[x_term],
        y=log_pvalue,
        hue=is_sig,
        palette={True: "purple", False: "lightgrey"},
        alpha=0.6,
    )

    plt.title(f"{title} - {target_variable.upper()} Effect")
    plt.xlabel("Effect Size (log2FC)")
    plt.ylabel("-log10(p-value)")
    plt.axhline(-np.log10(alpha), color="red", linestyle="--", alpha=0.5)

    safe_title = title.replace(" ", "_").replace("/", "-")
    fig_file = (
        figures_dir
        / f"{project}.{modality}.{regression_type}_volcano.{target_variable}.{safe_title}.png"
    )
    plt.savefig(fig_file, dpi=300)
    plt.close()
    logger.info("Saved volcano plot to %s", fig_file)


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figures_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Configure logging
    log_filename = (
        logs_dir
        / f"{args.modality}_{args.regression_type}_{args.target_variable}_post_regression.log"
    )
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logger.info("Command line: %s", " ".join(sys.argv))
    logger.info("Logging configured. Writing to %s", log_filename)

    project = args.project
    modality = args.modality.lower()
    regression_type = args.regression_type
    target_variable = args.target_variable

    # 1. Aggregation
    logger.info("Aggregating regression results...")

    # Pattern: {project}.{modality}.{cell_type}.{regression_type}.{target_variable}.csv
    pattern = f"{project}.{modality}.*.{regression_type}.{target_variable}.csv"
    result_files = list(results_dir.glob(pattern))

    if not result_files:
        logger.warning("No result files found matching %s in %s", pattern, results_dir)
        return

    dfs = []
    for file_path in result_files:
        try:
            df = read_csv(file_path)
            if "tissue" not in df.columns:
                continue
            dfs.append(df)
        except Exception as e:
            logger.error("Failed to read %s: %s", file_path, str(e))

    if not dfs:
        logger.error("No valid dataframes loaded.")
        return

    regression_results = concat(dfs, ignore_index=True)
    logger.info("Aggregated results shape: %s", regression_results.shape)

    if debug:
        print(regression_results.head())

    # 2. FDR Correction
    logger.info("Computing BH FDR...")
    regression_results["p-value"] = regression_results["p-value"].fillna(1)
    regression_results = compute_bh_fdr(regression_results, verbose=True)

    # 3. RLM Effect Size Filtering
    if regression_type == "rlm":
        logger.info(
            "Applying RLM effect size filter (min_effect=%.4f)...",
            args.min_rlm_effect,
        )
        mask = regression_results.coef.abs() < args.min_rlm_effect
        n_filtered = mask.sum()
        regression_results.loc[mask, "fdr_bh"] = 1.0
        logger.info("Set %d results to FDR=1.0 due to small effect size.", n_filtered)

    # Count significant
    n_sig = (regression_results["fdr_bh"] < 0.05).sum()
    logger.info("Total significant features (FDR < 0.05): %d", n_sig)

    # 4. Summary Counts
    total_counts = regression_results["tissue"].value_counts()
    sig_counts = regression_results.loc[regression_results["fdr_bh"] < 0.05][
        "tissue"
    ].value_counts()

    summary_df = pd.DataFrame({"Total": total_counts, "Significant": sig_counts})
    summary_df["Significant"] = summary_df["Significant"].fillna(0).astype(int)
    summary_df["Percentage"] = (
        summary_df["Significant"] / summary_df["Total"] * 100
    ).round(2)
    summary_df = summary_df.sort_values("Percentage", ascending=False)

    logger.info("Significant counts and percentages per cell type:\n%s", summary_df)

    # 5. Save Results
    results_file = (
        results_dir
        / f"{project}.all_celltypes.{modality}.{regression_type}.{target_variable}.csv"
    )
    results_fdr_file = (
        results_dir
        / f"{project}.all_celltypes.{modality}.{regression_type}_fdr.{target_variable}.csv"
    )

    logger.info("Saving full results to %s", results_file)
    regression_results.to_csv(results_file, index=False)

    logger.info("Saving significant results to %s", results_fdr_file)
    regression_results.loc[regression_results["fdr_bh"] < 0.05].to_csv(
        results_fdr_file, index=False
    )

    # 6. Volcano Plots
    logger.info("Generating volcano plots...")

    # Plot all results
    volcano_plot(
        regression_results,
        project,
        modality,
        regression_type,
        target_variable,
        figures_dir,
        title="All Cell Types",
    )

    # Plot per cell type
    if args.volcano_per_celltype:
        cell_types = regression_results["tissue"].unique()
        for ct in cell_types:
            ct_results = regression_results.loc[regression_results.tissue == ct]
            if ct_results.shape[0] > 0:
                volcano_plot(
                    ct_results,
                    project,
                    modality,
                    regression_type,
                    target_variable,
                    figures_dir,
                    title=ct,
                )


if __name__ == "__main__":
    main()
