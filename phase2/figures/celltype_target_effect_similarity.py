import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Visualize similarities between cell-types based on age-associated features."
    )
    parser.add_argument(
        "--project",
        type=str,
        default=DEFAULT_PROJECT,
        help="Project name used for file prefixes.",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default=DEFAULT_WRK_DIR,
        help="Base working directory.",
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
        help="Regression method to use.",
    )
    parser.add_argument(
        "--target-variable",
        type=str,
        default="age",
        help="The primary target variable column in covariates (default: 'age').",
    )
    parser.add_argument(
        "--effect-column",
        type=str,
        default="coef",
        choices=["coef", "z", "fc", "log2fc", "percentchange"],
        help="Effect column to use for Spearman correlation.",
    )
    parser.add_argument(
        "--feature-space",
        type=str,
        default="global",
        choices=["global", "pairwise"],
        help="Feature space for correlation: 'global' (union across all cell-types) or 'pairwise' (union per cell-type pair).",
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
    logs_dir = work_dir / "logs"
    figures_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Configure logging
    log_filename = (
        logs_dir / f"{args.modality}_{args.regression_type}_celltype_similarity.log"
    )
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Logging configured. Writing to {log_filename}")

    project = args.project
    modality = args.modality.lower()
    regression_type = args.regression_type
    effect_column = args.effect_column
    target_variable = args.target_variable
    feature_space = args.feature_space

    # Define file paths
    fdr_file = (
        results_dir
        / f"{project}.{modality}.all_celltypes.{regression_type}_fdr_filtered.{target_variable}.csv"
    )
    full_results_file = (
        results_dir / f"{project}.all_celltypes.{modality}.{regression_type}.{target_variable}.csv"
    )

    if not fdr_file.exists():
        logger.error(f"FDR filtered results file not found: {fdr_file}")
        return
    if not full_results_file.exists():
        logger.error(f"Full results file not found: {full_results_file}")
        return

    logger.info(f"Loading FDR filtered features from {fdr_file}")
    fdr_df = pd.read_csv(fdr_file)
    if "feature" not in fdr_df.columns:
        logger.error("Column 'feature' missing in FDR filtered results.")
        return

    sig_features = fdr_df["feature"].unique()
    logger.info(
        f"Found {len(sig_features)} unique significant features across all cell-types."
    )

    logger.info(f"Loading full results from {full_results_file}")
    full_df = pd.read_csv(full_results_file)

    if effect_column not in full_df.columns:
        logger.error(f"Effect column '{effect_column}' missing in full results.")
        return
    if "tissue" not in full_df.columns:
        logger.error("Column 'tissue' missing in full results.")
        return
    if "feature" not in full_df.columns:
        logger.error("Column 'feature' missing in full results.")
        return

    if feature_space == "global":
        logger.info(f"Filtering full results to {len(sig_features)} significant features.")
        filtered_df = full_df[full_df["feature"].isin(sig_features)].copy()

        # Pivot table
        logger.info(f"Pivoting table using effect column '{effect_column}'.")
        # Using drop_duplicates to handle any duplicated feature-tissue combinations
        pivot_df = filtered_df.drop_duplicates(subset=["feature", "tissue"]).pivot(
            index="feature", columns="tissue", values=effect_column
        )

        # Handle missing values
        missing_pct = pivot_df.isna().mean().mean() * 100
        if missing_pct > 0:
            logger.warning(
                f"Pivot table contains {missing_pct:.2f}% missing values. Filling with 0."
            )
            pivot_df = pivot_df.fillna(0)

        # Compute Spearman correlation
        logger.info("Computing Spearman correlation matrix (global feature space).")
        corr_matrix = pivot_df.corr(method="spearman")
    else:
        logger.info("Computing Spearman correlation matrix using pairwise feature space.")
        if "tissue" not in fdr_df.columns:
            logger.error(
                "Column 'tissue' missing in FDR filtered results; cannot compute pairwise feature space."
            )
            return

        full_pivot = full_df.drop_duplicates(subset=["feature", "tissue"]).pivot(
            index="feature", columns="tissue", values=effect_column
        )
        tissues = list(full_pivot.columns)
        sig_by_tissue = {
            t: set(fdr_df.loc[fdr_df["tissue"] == t, "feature"].unique())
            for t in tissues
        }

        corr_matrix = pd.DataFrame(index=tissues, columns=tissues, dtype=float)
        for i, t1 in enumerate(tissues):
            corr_matrix.loc[t1, t1] = 1.0
            for j in range(i + 1, len(tissues)):
                t2 = tissues[j]
                pair_union = list(sig_by_tissue[t1].union(sig_by_tissue[t2]))
                if len(pair_union) < 2:
                    val = 0.0
                else:
                    v1 = full_pivot[t1].reindex(pair_union).fillna(0)
                    v2 = full_pivot[t2].reindex(pair_union).fillna(0)
                    val = v1.corr(v2, method="spearman")
                    if pd.isna(val):
                        val = 0.0
                corr_matrix.loc[t1, t2] = val
                corr_matrix.loc[t2, t1] = val

    corr_matrix = corr_matrix.fillna(0)

    # Plot
    space_suffix = f".{feature_space}" if feature_space != "global" else ""
    fig_filename_png = (
        figures_dir
        / f"{project}.{modality}.{regression_type}.celltype_similarity.{target_variable}.{effect_column}{space_suffix}.png"
    )
    fig_filename_svg = (
        figures_dir
        / f"{project}.{modality}.{regression_type}.celltype_similarity.{target_variable}.{effect_column}{space_suffix}.svg"
    )
    logger.info("Generating clustered heatmap.")

    sns.set_theme(style="white")
    g = sns.clustermap(
        corr_matrix,
        cmap="vlag",
        annot=True,
        annot_kws={"size": 8},
        fmt=".2f",
        figsize=(10, 10),
        cbar_pos=(1.05, 0.2, 0.03, 0.6),
        vmin=-1,
        vmax=1,
        dendrogram_ratio=0.01,
    )
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    title = f"Cell-Type Similarity\nModality: {modality.upper()}, Target: {target_variable.upper()}, Effect: {effect_column}"
    if feature_space != "global":
        title += f", Space: {feature_space.capitalize()}"
    g.ax_heatmap.set_title(title, pad=20)
    g.ax_heatmap.set_xlabel("Cell Type")
    g.ax_heatmap.set_ylabel("Cell Type")

    g.savefig(fig_filename_png, dpi=300, bbox_inches="tight")
    g.savefig(fig_filename_svg, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved figures to {fig_filename_png} and {fig_filename_svg}")


if __name__ == "__main__":
    main()
