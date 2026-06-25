#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.stats.power import FTestPower


def compute_empirical_alpha(file_path):
    print(f"Reading WLS regression results from: {file_path} to compute alpha...")
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
    print(f"  -> Computed empirical alpha: {alpha:.4e}")
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


def process_dataset(label, n, results_path, var_threshold, alpha, target_power):
    print(f"\n" + "=" * 65)
    print(f"Dataset: {label} (N={n}, File: {results_path})")
    print(f"Empirical Alpha: {alpha:.4e}")
    print(
        f"Power-derived {int(target_power * 100)}% R2 Threshold: {var_threshold:.2f}%"
    )
    print("=" * 65)

    # 1. Load results file
    print(f"Reading WLS regression results from: {results_path}...")
    try:
        df = pd.read_csv(
            results_path,
            sep=",",
            usecols=["coef", "stderr", "p-value", "fdr_bh", "percentchange"],
        )
    except Exception as e:
        raise IOError(f"Failed to read/parse WLS result file '{results_path}': {e}")

    # 2. Calculate Variance Explained (R-squared)
    beta_sq = df["coef"] ** 2
    se_sq = df["stderr"] ** 2
    df["R_squared"] = beta_sq / (beta_sq + (se_sq * (n - 2)))
    df["Variance_Explained_Pct"] = df["R_squared"] * 100

    # Print top 5 by variance explained
    print("\nTop 5 Features by Variance Explained:")
    top_5 = df.sort_values("Variance_Explained_Pct", ascending=False).head()
    print(
        top_5[
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

    # Print bottom 5 significant ones by variance explained
    sig_df = df[df["fdr_bh"] <= 0.05]
    print("\nBottom 5 Significant Features (fdr_bh <= 0.05) by Variance Explained:")
    bottom_5 = sig_df.sort_values("Variance_Explained_Pct", ascending=True).head()
    print(
        bottom_5[
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

    # 3. Calculate threshold-based metrics (using dynamically solved threshold)
    total_sig = len(sig_df)
    sig_above_thresh = sig_df[sig_df["Variance_Explained_Pct"] >= var_threshold]
    count_above = len(sig_above_thresh)
    pct_above = (count_above / total_sig) * 100 if total_sig > 0 else 0.0

    # Find the empirical percent change at the R2 threshold
    df["abs_diff_from_thresh"] = np.abs(df["Variance_Explained_Pct"] - var_threshold)
    closest_row = df.loc[df["abs_diff_from_thresh"].idxmin()]
    empirical_pct_change_at_thresh = np.abs(closest_row["percentchange"])

    print("\n" + "-" * 50)
    print(
        f"Empirical Threshold from {int(target_power * 100)}% Power:"
    )
    print(f"  -> R-Squared: {var_threshold:.2f}%")
    print(f"  -> Approx Percent Change: {empirical_pct_change_at_thresh:.2f}%")
    print(f"Total Significant Features (fdr_bh <= 0.05): {total_sig:,}")
    print(f"Significant Features >= Threshold: {count_above:,} ({pct_above:.2f}%)")
    print("-" * 50)

    # 4. Generate rank-ordered S-curve scatter plot with downsampling to optimize plot speed and size
    print(f"Generating scatter plot of Variance Explained Pct for {label}...")
    df_sorted = df.sort_values("Variance_Explained_Pct", ascending=True).reset_index(
        drop=True
    )

    # Separate into significant and non-significant for plotting
    df_sig = df_sorted[df_sorted["fdr_bh"] <= 0.05]
    df_nonsig = df_sorted[df_sorted["fdr_bh"] > 0.05]

    # Systematic downsampling of non-significant points to keep plot fast and light (up to 10,000 points)
    max_nonsig_plot = 10000
    if len(df_nonsig) > max_nonsig_plot:
        indices = np.linspace(0, len(df_nonsig) - 1, max_nonsig_plot, dtype=int)
        df_nonsig_sampled = df_nonsig.iloc[indices]
    else:
        df_nonsig_sampled = df_nonsig

    x_nonsig = df_nonsig_sampled.index
    y_nonsig = df_nonsig_sampled["Variance_Explained_Pct"]

    # Systematic downsampling of significant points if they are too dense (up to 15,000 points)
    max_sig_plot = 15000
    if len(df_sig) > max_sig_plot:
        indices = np.linspace(0, len(df_sig) - 1, max_sig_plot, dtype=int)
        df_sig_sampled = df_sig.iloc[indices]
    else:
        df_sig_sampled = df_sig

    x_sig = df_sig_sampled.index
    y_sig = df_sig_sampled["Variance_Explained_Pct"]

    # Object-oriented matplotlib API to prevent memory leaks and ensure clean figure handles
    fig, ax = plt.subplots(figsize=(10, 6), dpi=300)

    # Plot non-significant points
    ax.scatter(
        x_nonsig,
        y_nonsig,
        color="#3b6a94",
        s=2,
        alpha=0.5,
        label="Features",
        edgecolors="none",
    )

    # Highlight significant features
    if len(x_sig) > 0:
        ax.scatter(
            x_sig,
            y_sig,
            color="#b05c60",
            s=3,
            alpha=0.7,
            label="Significant Features (fdr_bh <= 0.05)",
            edgecolors="none",
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
        f"Distribution of Variance Explained ($R^2$) across Age-Associated Features ({label}, N={n})",
        fontsize=14,
        fontweight="bold",
        pad=15,
    )
    ax.set_xlabel("Rank (sorted by Variance Explained)", fontsize=12, labelpad=10)
    ax.set_ylabel("Variance Explained ($R^2$) (%)", fontsize=12, labelpad=10)
    ax.grid(True, linestyle=":", alpha=0.6)
    ax.legend(loc="upper left", frameon=True, facecolor="white", edgecolor="none")

    plt.tight_layout()
    dist_fig_path = f"/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/figures/{label}_variance_explained_distribution.png"

    # Ensure output directory exists before saving
    os.makedirs(os.path.dirname(dist_fig_path), exist_ok=True)
    fig.savefig(dist_fig_path, dpi=300)
    plt.close(fig)
    print(f"Visualization successfully saved to: {dist_fig_path}")
    print("=" * 65 + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Unified tool for calculating statistical power for age-based WLS regressions and comparing with empirical variance explained."
    )
    parser.add_argument(
        "--labels",
        default="RNA,ATAC",
        help="Comma-separated dataset labels (default: 'RNA')",
    )
    parser.add_argument(
        "--sizes",
        default="35,35",
        help="Comma-separated dataset sample sizes (default: '100')",
    )
    parser.add_argument(
        "--results",
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.all_celltypes.rna.wls_fdr.age.csv,/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.all_celltypes.atac.wls_fdr.age.csv",
        help="Comma-separated paths to WLS regression result files (e.g. 'aging_phase2.all_celltypes.rna.wls_fdr.age.csv')",
    )
    parser.add_argument(
        "--target-power",
        type=float,
        default=0.80,
        help="Target power to solve for minimum detectable R-squared and empirical thresholding (default: 0.80)",
    )
    parser.add_argument(
        "--output",
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/figures/Age_WLS_Power_Curve.png",
        help="Path to save the generated theoretical power curves comparison plot (default: 'figures/Age_WLS_Power_Curve.png')",
    )
    args = parser.parse_args()

    # Parse comma-separated inputs
    labels = [x.strip() for x in args.labels.split(",")]
    sizes = [int(x.strip()) for x in args.sizes.split(",")]

    if not args.results:
        print("Please provide at least one result file using the --results argument.")
        return

    results_paths = [x.strip() for x in args.results.split(",")]
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
    print("=" * 65)
    print("Step 1: Computing empirical alphas from result sets...")
    print("=" * 65)
    for path in results_paths:
        alphas.append(compute_empirical_alpha(path))

    # Step 2: Solve for the theoretical minimum detectable R^2 for each dataset
    print("\n" + "=" * 65)
    print(
        f"Step 2: Solving for minimum detectable R^2 at {int(target_power * 100)}% power..."
    )
    print("=" * 65)
    var_thresholds = []
    for label, n, alpha in zip(labels, sizes, alphas):
        min_r2_pct = solve_min_detectable_r_squared(
            n, alpha, target_power, power_analysis
        )
        var_thresholds.append(min_r2_pct)
        print(f"  -> {label} minimum detectable R^2: {min_r2_pct:.2f}%")

    # Step 3: Run full empirical variance explained pipeline using the power-solved R^2 thresholds
    print("\n" + "=" * 65)
    print(
        "Step 3: Processing empirical dataset results using power-solved thresholds..."
    )
    print("=" * 65)
    for label, n, path, threshold, alpha in zip(
        labels, sizes, results_paths, var_thresholds, alphas
    ):
        process_dataset(label, n, path, threshold, alpha, target_power)

    # Step 4: Theoretical Power Curves comparison plot
    print("=" * 65)
    print("Step 4: Generating theoretical power comparison curves...")
    print("=" * 65)

    # Create array of R-squared values from 1% to 50%
    r_squared_range = np.linspace(0.01, 0.50, 100)
    effect_sizes_f_range = np.sqrt(r_squared_range / (1 - r_squared_range))

    # Object-oriented matplotlib API for final power curve plot
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    # Plot theoretical power curves
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

    # Formatting the plot
    ax.axhline(
        y=target_power,
        color="r",
        linestyle="--",
        label=f"{int(target_power * 100)}% Power Threshold",
    )
    ax.set_title(
        "Power to Detect Age Associations by Variance Explained ($R^2$)",
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

    # Save the comparison plot
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nTheoretical comparison power curve plot saved to: {args.output}")


if __name__ == "__main__":
    main()
