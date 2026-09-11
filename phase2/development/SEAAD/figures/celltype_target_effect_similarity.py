#!/usr/bin/env python3
"""
Visualize similarities of disease status effect sizes across cell types.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Configure logging
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize similarities of disease status effects across cell types."
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
        help="Regression method (e.g. wls, ols).",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


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
    log_filename = logs_dir / f"{args.modality}_{args.target_variable}_effect_similarity.log"
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

    # Load general filtered results
    general_fdr_file = (
        results_dir
        / f"{project}.{modality}.all_celltypes.{regression_type}_fdr_filtered.{target_variable}.csv"
    )

    logger.info("Reading filtered results from: %s", general_fdr_file)
    if not general_fdr_file.exists():
        logger.error("File not found: %s", general_fdr_file)
        sys.exit(1)
        
    results_df = pd.read_csv(general_fdr_file)
    
    if results_df.empty:
        logger.warning("No significant results found to compare. Exiting.")
        sys.exit(0)

    # We want to pivot the data so that rows are features, columns are cell types (tissue),
    # and values are the effect size (coef)
    logger.info("Pivoting effect sizes across cell types...")
    pivoted_df = results_df.pivot(index="feature", columns="tissue", values="coef")
    logger.info("Pivoted matrix shape: %s (features x cell types)", pivoted_df.shape)

    # Compute correlation matrix between cell types
    # Pearson correlation is typical for comparing continuous effect size values
    corr_matrix = pivoted_df.corr(method="pearson").fillna(0)

    # Generate clustered heatmap
    try:
        n_cell_types = corr_matrix.shape[0]
        fig_size = max(8, n_cell_types * 0.4)
        
        plt.figure(figsize=(fig_size, fig_size))
        g = sns.clustermap(
            corr_matrix,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            center=0,
            figsize=(fig_size, fig_size),
            xticklabels=True,
            yticklabels=True,
            linewidths=0.5,
            annot=True,
            fmt=".2f",
            cbar_kws={"label": "Pearson Correlation"},
        )
        
        g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=9, rotation=90)
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=9, rotation=0)
        
        # Add titles
        plt.suptitle(
            f"Similarity of {target_variable.upper()} Effect Sizes across Cell Types ({modality.upper()})\n({regression_type.upper()} model, Pearson correlation)",
            fontsize=12,
            fontweight="bold",
            y=1.02
        )
        
        heatmap_out_png = (
            figures_dir / f"{project}_{modality}_{target_variable}_effect_similarity_clustermap.png"
        )
        g.savefig(heatmap_out_png, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info("Saved clustered similarity heatmap to %s", heatmap_out_png)

    except Exception as e:
        logger.error("Failed to generate similarity clustermap: %s", str(e))


if __name__ == "__main__":
    main()
