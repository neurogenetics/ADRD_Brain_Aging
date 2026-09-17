#!/usr/bin/env python3
"""
Compare Age analysis results with Dx (disease-status) analysis results.
Computes the intersection of significant features and generates a cross-correlation similarity heatmap.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from tabulate import tabulate

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare Age and Dx regression results across cell types and modalities."
    )
    parser.add_argument(
        "--age-project",
        type=str,
        default="aging_phase2",
        help="Project name used for the Age analysis file prefixes.",
    )
    parser.add_argument(
        "--age-work-dir",
        type=str,
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2",
        help="Base working directory for the Age analysis.",
    )
    parser.add_argument(
        "--dx-project",
        type=str,
        default="seaad_ec_multiome",
        help="Project name used for the Dx analysis file prefixes.",
    )
    parser.add_argument(
        "--dx-work-dir",
        type=str,
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad",
        help="Base working directory for the Dx analysis.",
    )
    parser.add_argument(
        "--modality",
        type=str,
        required=True,
        choices=["rna", "atac"],
        help="Data modality (rna or atac).",
    )
    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
        help="Regression method to use (e.g. wls, ols).",
    )
    parser.add_argument(
        "--effect-column",
        type=str,
        default="coef",
        choices=["coef", "z", "fc", "log2fc", "percentchange"],
        help="Effect column to use for Spearman correlation.",
    )
    parser.add_argument(
        "--cell-type-map",
        type=str,
        default=None,
        help="Comma-separated mapping of Age cell-type names to Dx names, format: 'AgeName:DxName,AgeName2:DxName2'",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def parse_cell_type_map(map_str):
    if not map_str:
        return {}
    mapping = {}
    for item in map_str.split(","):
        if ":" in item:
            k, v = item.split(":", 1)
            mapping[k.strip()] = v.strip()
    return mapping


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    age_work_dir = Path(args.age_work_dir)
    age_results_dir = age_work_dir / "results"

    dx_work_dir = Path(args.dx_work_dir)
    dx_results_dir = dx_work_dir / "results"
    dx_figures_dir = dx_work_dir / "figures"
    dx_logs_dir = dx_work_dir / "logs"

    dx_figures_dir.mkdir(parents=True, exist_ok=True)
    dx_logs_dir.mkdir(parents=True, exist_ok=True)

    # Configure logging to file and stdout
    log_filename = dx_logs_dir / f"{args.modality}_compare_age_dx_effects.log"
    file_handler = logging.FileHandler(log_filename, mode="w")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting Age vs Dx Comparison ===")
    logger.info("Logging configured. Writing to %s", log_filename)

    modality = args.modality.lower()
    regression_type = args.regression_type
    effect_column = args.effect_column

    # Define input file paths
    # Age inputs
    age_fdr_file = (
        age_results_dir
        / f"{args.age_project}.{modality}.all_celltypes.{regression_type}_fdr_filtered.age.csv"
    )
    age_full_file = (
        age_results_dir / f"{args.age_project}.all_celltypes.{modality}.{regression_type}.age.csv"
    )

    # Dx inputs
    dx_fdr_file = (
        dx_results_dir
        / f"{args.dx_project}.{modality}.all_celltypes.{regression_type}_fdr_filtered.dx.csv"
    )
    dx_full_file = (
        dx_results_dir / f"{args.dx_project}.all_celltypes.{modality}.{regression_type}.dx.csv"
    )

    # Check file existence
    files_missing = False
    for label, path in [
        ("Age FDR", age_fdr_file),
        ("Age Full", age_full_file),
        ("Dx FDR", dx_fdr_file),
        ("Dx Full", dx_full_file),
    ]:
        if not path.exists():
            logger.error("%s file not found at: %s", label, path)
            files_missing = True
    if files_missing:
        sys.exit(1)

    # 1. Load data
    logger.info("Loading results...")
    age_fdr = pd.read_csv(age_fdr_file)
    age_full = pd.read_csv(age_full_file)
    dx_fdr = pd.read_csv(dx_fdr_file)
    dx_full = pd.read_csv(dx_full_file)

    # Apply cell-type name mapping to Age data if provided
    cell_type_map = parse_cell_type_map(args.cell_type_map)
    if cell_type_map:
        logger.info("Applying cell-type mapping to Age data: %s", cell_type_map)
        if "tissue" in age_fdr.columns:
            age_fdr["tissue"] = age_fdr["tissue"].replace(cell_type_map)
        if "tissue" in age_full.columns:
            age_full["tissue"] = age_full["tissue"].replace(cell_type_map)

    # 2. Compute significant feature intersections per cell type (tissue)
    logger.info("Computing significant feature overlaps per cell type...")
    
    # Common cell types (intersection of tissue names in both datasets)
    age_tissues = set(age_fdr["tissue"].unique()) if "tissue" in age_fdr.columns else set()
    dx_tissues = set(dx_fdr["tissue"].unique()) if "tissue" in dx_fdr.columns else set()
    common_tissues = sorted(list(age_tissues.intersection(dx_tissues)))

    overlap_results = []
    for tissue in common_tissues:
        age_sig = set(age_fdr.loc[age_fdr["tissue"] == tissue, "feature"].unique())
        dx_sig = set(dx_fdr.loc[dx_fdr["tissue"] == tissue, "feature"].unique())
        shared_sig = age_sig.intersection(dx_sig)
        
        union_len = len(age_sig.union(dx_sig))
        jaccard = (len(shared_sig) / union_len) if union_len > 0 else 0.0
        
        dx_len = len(dx_sig)
        pct_dx_age_associated = (len(shared_sig) / dx_len * 100) if dx_len > 0 else 0.0

        overlap_results.append({
            "tissue": tissue,
            "age_significant_count": len(age_sig),
            "dx_significant_count": len(dx_sig),
            "shared_significant_count": len(shared_sig),
            "jaccard_similarity": round(jaccard, 4),
            "percent_dx_features_age_associated": round(pct_dx_age_associated, 2),
        })

    overlap_df = pd.DataFrame(overlap_results)
    
    # Save overlap table to CSV
    overlap_csv = dx_results_dir / f"{args.dx_project}_{modality}_{regression_type}_age_dx_sig_overlap.csv"
    overlap_df.to_csv(overlap_csv, index=False)
    logger.info("Saved significant overlap metrics to %s", overlap_csv)
    
    print("\n=== Significant Feature Overlap Summary ===")
    print(tabulate(overlap_df, headers="keys", tablefmt="psql", showindex=False))

    # 3. Compute cross-correlation of effect sizes (clustermap)
    # Get the union of significant features across both age and dx to use as background
    age_sig_features = set(age_fdr["feature"].unique())
    dx_sig_features = set(dx_fdr["feature"].unique())
    sig_features = sorted(list(age_sig_features.union(dx_sig_features)))

    logger.info("Union of significant features across both: %d", len(sig_features))
    
    # Filter full results to these features
    age_filtered = age_full[age_full["feature"].isin(sig_features)].copy()
    dx_filtered = dx_full[dx_full["feature"].isin(sig_features)].copy()

    # Pivot tables
    logger.info("Pivoting tables...")
    age_pivot = age_filtered.drop_duplicates(subset=["feature", "tissue"]).pivot(
        index="feature", columns="tissue", values=effect_column
    )
    dx_pivot = dx_filtered.drop_duplicates(subset=["feature", "tissue"]).pivot(
        index="feature", columns="tissue", values=effect_column
    )

    # Align indexes on common features
    common_idx = age_pivot.index.intersection(dx_pivot.index)
    if len(common_idx) == 0:
        logger.warning("No overlapping features in full tables. Heatmap generation skipped.")
        sys.exit(0)

    logger.info("Aligning %d features for correlation matrix...", len(common_idx))
    age_aligned = age_pivot.loc[common_idx]
    dx_aligned = dx_pivot.loc[common_idx]

    # Combine columns under keys to compute full correlation
    combined = pd.concat([age_aligned, dx_aligned], axis=1, keys=["Age", "Dx"])
    
    # Handle missing values
    missing_pct = combined.isna().mean().mean() * 100
    if missing_pct > 0:
        logger.warning("Combined matrix contains %.2f%% missing values. Filling with 0.", missing_pct)
        combined = combined.fillna(0)

    # Compute correlation
    logger.info("Computing Spearman cross-correlation...")
    corr = combined.corr(method="spearman")
    
    # Extract asymmetric sub-matrix: Rows = Dx cell types, Columns = Age cell types
    cross_corr = corr.loc["Dx", "Age"].fillna(0)

    # Plot
    fig_filename_png = (
        dx_figures_dir
        / f"{args.dx_project}.{modality}.{regression_type}_age_dx_similarity.{effect_column}.png"
    )
    fig_filename_svg = (
        dx_figures_dir
        / f"{args.dx_project}.{modality}.{regression_type}_age_dx_similarity.{effect_column}.svg"
    )

    logger.info("Generating asymmetric clustered heatmap...")
    try:
        # Clustermap handles asymmetric matrices beautifully
        plt.figure(figsize=(10, 8))
        sns.set_theme(style="white")
        
        g = sns.clustermap(
            cross_corr,
            cmap="vlag",
            annot=True,
            annot_kws={"size": 8},
            fmt=".2f",
            figsize=(10, 8),
            vmin=-1,
            vmax=1,
            cbar_pos=(1.05, 0.2, 0.03, 0.6),
            dendrogram_ratio=0.01,
        )
        g.ax_row_dendrogram.set_visible(False)
        g.ax_col_dendrogram.set_visible(False)
        
        g.ax_heatmap.set_title(
            f"Age vs Dx Cell-Type Similarity\nModality: {modality.upper()}, Effect: {effect_column}",
            pad=20,
            fontweight="bold",
        )
        g.ax_heatmap.set_xlabel("Age Cell Types (aging_phase2)")
        g.ax_heatmap.set_ylabel("Dx Cell Types (seaad_ec_multiome)")

        g.savefig(fig_filename_png, dpi=300, bbox_inches="tight")
        g.savefig(fig_filename_svg, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("Saved comparison heatmap to %s and %s", fig_filename_png, fig_filename_svg)
    except Exception as e:
        logger.error("Failed to generate similarity clustermap: %s", str(e))

    logger.info("=== Age vs Dx Comparison Finished Successfully ===")


if __name__ == "__main__":
    main()
