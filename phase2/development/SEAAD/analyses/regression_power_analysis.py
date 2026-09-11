#!/usr/bin/env python3
"""
Unified tool for calculating statistical power for disease-based regressions and comparing with empirical variance explained.
"""

import os
import argparse
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.power import FTestPower

# Setup module-level logger
logger = logging.getLogger(__name__)


def compute_empirical_alpha(file_path):
    logger.info("Reading WLS regression results from: %s to compute alpha...", file_path)
    try:
        # Read only the necessary columns to save memory and keep execution fast
        df = pd.read_csv(file_path, sep=",", usecols=["p-value", "fdr_bh"])
    except Exception as e:
        raise IOError(f"Failed to read/parse WLS result file '{file_path}': {e}")

    sig_df = df[df["fdr_bh"] <= 0.05]
    if sig_df.empty:
        raise ValueError(
            f"No significant associations (fdr_bh <= 0.05) found in {file_path} to compute alpha."
        )

    alpha = sig_df["p-value"].max()
    logger.info("  -> Computed empirical alpha: %.4e", alpha)
    return alpha


def solve_min_detectable_r_squared(n, alpha, target_power, power_analysis):
    effect_size_f = power_analysis.solve_power(
        effect_size=None,
        df_num=n - 2,
        df_denom=1,
        alpha=alpha,
        power=target_power,
        ncc=1,
    )
    f_squared = effect_size_f**2
    r_squared = f_squared / (1 + f_squared)

    return r_squared * 100


def plot_variance_distribution(df, label, n, var_threshold, target_power, output_dir, target_variable):
    """Generates and saves a rank-ordered S-curve scatter plot of variance explained."""
    logger.info("Generating scatter plot of Variance Explained Pct for %s...", label)
    df_sorted = df.sort_values("Variance_Explained_Pct", ascending=True).reset_index(
        drop=True
    )

    # Separate into significant and non-significant for plotting
    df_sig = df_sorted[df_sorted["fdr_bh"] <= 0.05]
    df_nonsig = df_sorted[df_sorted["fdr_bh"] > 0.05]

    # Downsample non-significant points to keep plot light
    max_nonsig_plot = 10000
    if len(df_nonsig) > max_nonsig_plot:
        indices = np.linspace(0, len(df_nonsig) - 1, max_nonsig_plot, dtype=int)
        df_nonsig_sampled = df_nonsig.iloc[indices]
    else:
        df_nonsig_sampled = df_nonsig

    x_nonsig = df_nonsig_sampled.index
    y_nonsig = df_nonsig_sampled["Variance_Explained_Pct"]

    # Downsample significant points
    max_sig_plot = 15000
    if len(df_sig) > max_sig_plot:
        indices = np.linspace(0, len(df_sig) - 1, max_sig_plot, dtype=int)
        df_sig_sampled = df_sig.iloc[indices]
    else:
        df_sig_sampled = df_sig

    x_sig = df_sig_sampled.index
    y_sig = df_sig_sampled["Variance_Explained_Pct"]

    # Matplotlib OO API
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # Plot non-significant points
    ax.scatter(
        x_nonsig,
        y_nonsig,
        color="#3b6a94",
        s=1.5,
        alpha=0.4,
        label="Non-significant Features",
        edgecolors="none",
        zorder=1,
    )

    # Highlight significant features
    if len(x_sig) > 0:
        ax.scatter(
            x_sig,
            y_sig,
            color="#b05c60",
            s=0.5,
            alpha=0.6,
            label="Significant Features (fdr_bh <= 0.05)",
            edgecolors="none",
            zorder=2,
        )

    # Add horizontal threshold line
    ax.axhline(
        y=var_threshold,
        color="#d62728",
        linestyle="--",
        linewidth=1.5,
        label=f"{int(target_power * 100)}% Power Threshold ({var_threshold:.2f}%)",
    )

    ax.set_title(
        f"Distribution of Variance Explained ($R^2$) across {target_variable.upper()}-Associated Features ({label}, N={n})",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    ax.set_xlabel("Rank (sorted by Variance Explained)", fontsize=12, labelpad=10)
    ax.set_ylabel("Variance Explained ($R^2$) (%)", fontsize=12, labelpad=10)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none")

    plt.tight_layout()
    dist_fig_path = os.path.join(
        output_dir, f"{label}_variance_explained_distribution.{target_variable}.png"
    )

    fig.savefig(dist_fig_path, dpi=300)
    plt.close(fig)
    logger.info("Visualization successfully saved to: %s", dist_fig_path)


def process_dataset(
    label, n, results_path, var_threshold, alpha, target_power, output_dir, target_variable
):
    logger.info("\n" + "=" * 65)
    logger.info("Dataset: %s (N=%d, File: %s)", label, n, results_path)
    logger.info("Empirical Alpha: %.4e", alpha)
    logger.info(
        "Power-derived %d%% R2 Threshold: %.2f%%", int(target_power * 100), var_threshold
    )
    logger.info("=" * 65)

    # 1. Load results file
    logger.info("Reading regression results from: %s...", results_path)
    try:
        df = pd.read_csv(
            results_path,
            sep=",",
            usecols=["tissue", "coef", "stderr", "p-value", "fdr_bh", "percentchange"],
        )
    except Exception as e:
        raise IOError(f"Failed to read/parse result file '{results_path}': {e}")

    # 2. Calculate Variance Explained (R-squared)
    beta_sq = df["coef"] ** 2
    se_sq = df["stderr"] ** 2
    df["R_squared"] = beta_sq / (beta_sq + (se_sq * (n - 2)))
    df["Variance_Explained_Pct"] = df["R_squared"] * 100

    # Print top 5
    logger.info("\nTop 5 Features by Variance Explained:")
    top_5 = df.sort_values("Variance_Explained_Pct", ascending=False).head()
    logger.info(
        "\n"
        + top_5[
            [
                "coef",
                "percentchange",
                "stderr",
                "p-value",
                "fdr_bh",
                "Variance_Explained_Pct",
            ]
        ].to_string(index=False)
    )

    # Print bottom 5 significant ones
    sig_df = df[df["fdr_bh"] <= 0.05]
    logger.info(
        "\nBottom 5 Significant Features (fdr_bh <= 0.05) by Variance Explained:"
    )
    bottom_5 = sig_df.sort_values("Variance_Explained_Pct", ascending=True).head()
    logger.info(
        "\n"
        + bottom_5[
            [
                "coef",
                "percentchange",
                "stderr",
                "p-value",
                "fdr_bh",
                "Variance_Explained_Pct",
            ]
        ].to_string(index=False)
    )

    # 3. Calculate threshold-based metrics
    total_sig = len(sig_df)
    sig_above_thresh = sig_df[sig_df["Variance_Explained_Pct"] >= var_threshold]
    count_above = len(sig_above_thresh)
    pct_above = (count_above / total_sig) * 100 if total_sig > 0 else 0.0

    # Find empirical percent change at the R2 threshold
    df["abs_diff_from_thresh"] = np.abs(df["Variance_Explained_Pct"] - var_threshold)
    closest_row = df.loc[df["abs_diff_from_thresh"].idxmin()]
    empirical_pct_change_at_thresh = np.abs(closest_row["percentchange"])

    logger.info("\n" + "-" * 50)
    logger.info("Empirical Threshold from %d%% Power:", int(target_power * 100))
    logger.info("  -> R-Squared: %.2f%%", var_threshold)
    logger.info("  -> Approx Percent Change: %.2f%%", empirical_pct_change_at_thresh)
    logger.info("Total Significant Features (fdr_bh <= 0.05): %d", total_sig)
    logger.info(
        "Significant Features >= Threshold: %d (%.2f%%)", count_above, pct_above
    )

    logger.info("\nSignificant Features >= Threshold per Cell-Type (tissue):")
    if "tissue" in sig_df.columns:
        tissue_counts = sig_df["tissue"].value_counts()
        tissue_above_counts = sig_above_thresh["tissue"].value_counts()
        for tissue, total in tissue_counts.items():
            above = tissue_above_counts.get(tissue, 0)
            pct = (above / total) * 100 if total > 0 else 0.0
            logger.info("  -> %s: %d / %d (%.2f%%)", tissue, above, total, pct)
    else:
        logger.info("  -> 'tissue' column not found in data.")

    logger.info("-" * 50)

    # 4. Generate S-curve scatter plot
    plot_variance_distribution(df, label, n, var_threshold, target_power, output_dir, target_variable)


def main():
    parser = argparse.ArgumentParser(
        description="Unified tool for calculating statistical power for disease-based regressions and comparing with empirical variance."
    )
    parser.add_argument(
        "--project",
        default="seaad_ec_multiome",
        help="Project name used for file prefixes.",
    )
    parser.add_argument(
        "--work-dir",
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad",
        help="Base working directory.",
    )
    parser.add_argument(
        "--target-variable",
        default="dx",
        help="The primary disease target variable column in covariates (default: 'dx').",
    )
    parser.add_argument(
        "--regression-type",
        default="ols",
        help="Regression type (e.g. ols, wls) to examine in results directory.",
    )
    parser.add_argument(
        "--labels",
        default="RNA,ATAC",
        help="Comma-separated dataset labels (default: 'RNA,ATAC')",
    )
    parser.add_argument(
        "--sizes",
        default="35,35",
        help="Comma-separated dataset sample sizes (default: '35,35')",
    )
    parser.add_argument(
        "--results",
        default=None,
        help="Comma-separated paths to WLS regression result files. If omitted, they are resolved dynamically from output directories.",
    )
    parser.add_argument(
        "--target-power",
        type=float,
        default=0.80,
        help="Target power to solve for minimum detectable R-squared (default: 0.80)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Path to save the generated theoretical power curves. Defaults to figures/[target_variable]_WLS_Power_Curve.png",
    )
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figs_dir = work_dir / "figures"
    figs_dir.mkdir(parents=True, exist_ok=True)

    output_path = args.output if args.output else str(figs_dir / f"{args.target_variable}_WLS_Power_Curve.png")
    output_dir = os.path.dirname(output_path) or "."
    
    # Set up File Logging next to plot outputs
    log_file_path = os.path.join(output_dir, f"{args.target_variable}_regression_power_analysis.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file_path, mode="w"),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )

    logger.info("=" * 65)
    logger.info("Executing disease regression_power_analysis.py:")
    logger.info("  --project:         %s", args.project)
    logger.info("  --target-variable: %s", args.target_variable)
    logger.info("  --labels:          %s", args.labels)
    logger.info("  --sizes:           %s", args.sizes)
    logger.info("  --target-power:    %s", args.target_power)
    logger.info("  --output:          %s", output_path)
    logger.info("=" * 65 + "\n")

    # Resolve default results paths if none provided
    if args.results:
        results_paths = [x.strip() for x in args.results.split(",")]
    else:
        results_paths = [
            str(results_dir / f"{args.project}.all_celltypes.rna.{args.regression_type}_fdr.{args.target_variable}.csv"),
            str(results_dir / f"{args.project}.all_celltypes.atac.{args.regression_type}_fdr.{args.target_variable}.csv"),
        ]

    labels = [x.strip() for x in args.labels.split(",")]
    sizes = [int(x.strip()) for x in args.sizes.split(",")]
    target_power = args.target_power

    if not (len(labels) == len(sizes) == len(results_paths)):
        raise ValueError(
            f"The number of labels ({len(labels)}), sizes ({len(sizes)}), "
            f"and result files ({len(results_paths)}) must be equal."
        )

    # Initialize the power analysis object
    power_analysis = FTestPower()

    # Step 1: Compute empirical alphas dynamically from result files
    alphas = []
    logger.info("=" * 65)
    logger.info("Step 1: Computing empirical alphas from result sets...")
    logger.info("=" * 65)
    for path in results_paths:
        alphas.append(compute_empirical_alpha(path))

    # Step 2: Solve for theoretical minimum detectable R^2
    logger.info("\n" + "=" * 65)
    logger.info(
        f"Step 2: Solving for minimum detectable R^2 at {int(target_power * 100)}% power..."
    )
    logger.info("=" * 65)
    var_thresholds = []
    for label, n, alpha in zip(labels, sizes, alphas):
        min_r2_pct = solve_min_detectable_r_squared(
            n, alpha, target_power, power_analysis
        )
        var_thresholds.append(min_r2_pct)
        logger.info(f"  -> {label} minimum detectable R^2: {min_r2_pct:.2f}%")

    # Step 3: Process empirical variance explained
    logger.info("\n" + "=" * 65)
    logger.info(
        "Step 3: Processing empirical dataset results using power-solved thresholds..."
    )
    logger.info("=" * 65)
    for label, n, path, threshold, alpha in zip(
        labels, sizes, results_paths, var_thresholds, alphas
    ):
        process_dataset(label, n, path, threshold, alpha, target_power, output_dir, args.target_variable)

    # Step 4: Theoretical Power Curves comparison plot
    logger.info("=" * 65)
    logger.info("Step 4: Generating theoretical power comparison curves...")
    logger.info("=" * 65)

    r_squared_range = np.linspace(0.01, 0.50, 100)
    effect_sizes_f_range = np.sqrt(r_squared_range / (1 - r_squared_range))

    # Matplotlib OO API
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Plot curves
    for label, n, alpha in zip(labels, sizes, alphas):
        powers = power_analysis.power(
            effect_size=effect_sizes_f_range,
            df_num=n - 2,
            df_denom=1,
            alpha=alpha,
            ncc=1,
        )
        ax.plot(
            r_squared_range * 100,
            powers,
            label=f"{label} (n={n}, $\\alpha$={alpha:.1e})",
            linewidth=2,
        )

    # Formatting plot
    ax.axhline(
        y=target_power,
        color="r",
        linestyle="--",
        label=f"{int(target_power * 100)}% Power Threshold",
    )
    ax.set_title(
        f"Power to Detect {args.target_variable.upper()} Associations by Variance Explained ($R^2$)",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    ax.set_xlabel("Variance Explained ($R^2$) (%)", fontsize=12, labelpad=10)
    ax.set_ylabel("Statistical Power", fontsize=12, labelpad=10)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(0, max(r_squared_range * 100))
    ax.legend(loc="lower right", fontsize=11)
    ax.grid(True, linestyle=":", alpha=0.7)

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logger.info("\nTheoretical comparison power curve plot saved to: %s", output_path)


if __name__ == "__main__":
    main()
