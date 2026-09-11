#!/usr/bin/env python3
"""
Filter regression results based on comparison between General (LM) and Robust (RLM) models.
"""

import sys
import logging
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats

# Configure logging
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter regression results based on comparison between General (LM) and Robust (RLM) models."
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
        "--general-type",
        type=str,
        default="wls",
        help="General regression method (e.g., wls, ols).",
    )
    parser.add_argument(
        "--robust-type",
        type=str,
        default="vwrlm",
        help="Robust regression method (e.g., vwrlm, rlm).",
    )
    parser.add_argument(
        "--max-z",
        type=float,
        default=3.0,
        help="Maximum Z-score for effect size difference outlier filtering.",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def create_pair_id(df: pd.DataFrame) -> pd.Series:
    """Create a unique identifier for feature-tissue pairs."""
    if "endo_feature" in df.columns and "exog_feature" in df.columns:
        return df["endo_feature"] + "_" + df["exog_feature"] + "_" + df["tissue"]
    elif "feature" in df.columns and "tissue" in df.columns:
        return df["feature"] + "_" + df["tissue"]
    else:
        raise ValueError(
            "DataFrame must contain ('feature', 'tissue') or ('endo_feature', 'exog_feature', 'tissue') columns."
        )


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figures_dir = work_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    logs_dir = work_dir / "logs"

    # Configure logging
    log_filename = f"{logs_dir}/{args.modality}_{args.target_variable}_filter_regression_differences.log"
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
    general_type = args.general_type
    robust_type = args.robust_type
    target_variable = args.target_variable

    # Define file paths using dynamic target variable
    general_fdr_file = (
        results_dir
        / f"{project}.all_celltypes.{modality}.{general_type}_fdr.{target_variable}.csv"
    )
    robust_file = (
        results_dir
        / f"{project}.all_celltypes.{modality}.{robust_type}.{target_variable}.csv"
    )
    results_file = (
        results_dir
        / f"{project}.{modality}.all_celltypes.{general_type}_fdr_filtered.{target_variable}.csv"
    )

    logger.info("Reading general results (LM FDR filtered) from: %s", general_fdr_file)
    if not general_fdr_file.exists():
        logger.error("File not found: %s", general_fdr_file)
        return
    general_results = pd.read_csv(general_fdr_file)
    logger.info("General results shape: %s", general_results.shape)

    logger.info("Reading robust results (RLM Full) from: %s", robust_file)
    if not robust_file.exists():
        logger.error("File not found: %s", robust_file)
        return
    robust_results = pd.read_csv(robust_file)
    logger.info("Robust results shape: %s", robust_results.shape)

    if debug:
        print("General sample:")
        print(general_results.head())
        print("Robust sample:")
        print(robust_results.head())

    # 1. Filter Robust Results (Nominal p <= 0.05)
    logger.info("Filtering robust results for p-value <= 0.05...")
    robust_results = robust_results.loc[robust_results["p-value"] <= 0.05].copy()
    logger.info(
        "Robust results shape after p-value filtering: %s", robust_results.shape
    )

    # 2. Create Pair IDs
    logger.info("Creating pair IDs...")
    general_results["pair"] = create_pair_id(general_results)
    robust_results["pair"] = create_pair_id(robust_results)

    # 3. Find Intersection
    pair_intersect = set(general_results.pair) & set(robust_results.pair)
    pct_intersect = (
        (len(pair_intersect) / general_results.shape[0]) * 100
        if general_results.shape[0] > 0
        else 0
    )
    logger.info(
        "Intersection: %d features found in both (%.2f%%)",
        len(pair_intersect),
        pct_intersect,
    )

    # Filter general results to intersection
    filtered_results = general_results.loc[
        general_results.pair.isin(pair_intersect)
    ].copy()

    # 4. Check Direction Consistency
    logger.info("Checking effect direction consistency...")
    merged = filtered_results.merge(
        robust_results[["pair", "coef"]],
        on="pair",
        how="left",
        suffixes=(f"_{general_type}", f"_{robust_type}"),
    )

    coef_gen_col = f"coef_{general_type}"
    coef_rob_col = f"coef_{robust_type}"

    # Identify retained pairs (same direction)
    consistent_mask = (merged[coef_gen_col] * merged[coef_rob_col]) >= 0
    kept_pairs = merged.loc[consistent_mask, "pair"]

    logger.info(
        "Features with consistent direction: %d out of %d",
        len(kept_pairs),
        len(merged),
    )

    # Filter
    filtered_results = filtered_results.loc[
        filtered_results.pair.isin(kept_pairs)
    ].copy()
    kept_merged = merged.loc[consistent_mask].copy()

    # 5. Outlier Filtering (Z-score of effect delta)
    logger.info("Filtering effect size outliers (Max Z = %.2f)...", args.max_z)
    kept_merged["effect_delta"] = kept_merged[coef_gen_col] - kept_merged[coef_rob_col]

    # Avoid z-score errors if we only have 1 or 2 entries in testing
    if len(kept_merged) > 1:
        kept_merged["effect_delta_z"] = stats.zscore(kept_merged["effect_delta"])
    else:
        kept_merged["effect_delta_z"] = 0.0

    # Plot distribution
    try:
        plt.figure(figsize=(10, 6))
        sns.histplot(kept_merged["effect_delta_z"], kde=True)
        plt.title(
            f"Distribution of Effect Delta Z-scores ({modality} - {target_variable.upper()})"
        )
        plt.xlabel(
            f"Z-score of (Coef {general_type.upper()} - Coef {robust_type.upper()})"
        )
        fig_path = (
            figures_dir
            / f"{project}.{modality}.effect_delta_z_distribution.{target_variable}.png"
        )
        plt.savefig(fig_path)
        plt.close()
        logger.info("Saved Z-score distribution plot to %s", fig_path)
    except Exception as e:
        logger.warning("Failed to generate Z-score distribution plot: %s", str(e))

    # Apply filter
    z_mask = kept_merged["effect_delta_z"].abs() < args.max_z
    final_pairs = kept_merged.loc[z_mask, "pair"]

    n_outliers = (~z_mask).sum()
    logger.info("Removed %d outliers based on Z-score.", n_outliers)

    final_results = filtered_results.loc[filtered_results.pair.isin(final_pairs)].copy()

    # 6. Save
    logger.info("Final filtered results shape: %s", final_results.shape)

    # Drop the 'pair' column before saving to match original structure
    final_results = final_results.drop(columns=["pair"])

    logger.info("Saving results to %s", results_file)
    final_results.to_csv(results_file, index=False)

    # Summary of counts per cell type
    total_counts = general_results["tissue"].value_counts()
    sig_counts = final_results["tissue"].value_counts()

    summary_df = pd.DataFrame({"Total": total_counts, "Retained": sig_counts})
    summary_df["Retained"] = summary_df["Retained"].fillna(0).astype(int)
    summary_df["Percentage"] = (
        summary_df["Retained"] / summary_df["Total"] * 100
    ).round(2)
    summary_df = summary_df.sort_values("Percentage", ascending=False)

    logger.info(
        "Summary of filtered results per cell type:\n%s", summary_df.sort_index()
    )


if __name__ == "__main__":
    main()
