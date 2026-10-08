import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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
    parser.add_argument(
        "--plot-type",
        type=str,
        default="heatmap",
        choices=["heatmap", "dotplot"],
        help="Visualization type: 'heatmap' (default) or 'dotplot' (circle diameter reflects shared feature count).",
    )
    parser.add_argument(
        "--include-diagonal",
        action="store_true",
        help="Include main diagonal (self-comparison) in dotplot visualization (default: False).",
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
    plot_type = args.plot_type
    include_diagonal = args.include_diagonal

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

    # Compute intersection count matrix for dotplot or reporting
    tissues = list(corr_matrix.columns)
    if "tissue" in fdr_df.columns:
        sig_by_tissue_map = {
            t: set(fdr_df.loc[fdr_df["tissue"] == t, "feature"].unique())
            for t in tissues
        }
    else:
        sig_by_tissue_map = {t: set() for t in tissues}

    count_matrix = pd.DataFrame(index=tissues, columns=tissues, dtype=int)
    for i, t1 in enumerate(tissues):
        if include_diagonal:
            count_matrix.loc[t1, t1] = len(sig_by_tissue_map.get(t1, set()))
        else:
            count_matrix.loc[t1, t1] = 0
        for j in range(i + 1, len(tissues)):
            t2 = tissues[j]
            intersection_count = len(
                sig_by_tissue_map.get(t1, set()).intersection(sig_by_tissue_map.get(t2, set()))
            )
            count_matrix.loc[t1, t2] = intersection_count
            count_matrix.loc[t2, t1] = intersection_count

    # Plot
    space_suffix = f".{feature_space}" if feature_space != "global" else ""
    plot_suffix = f".{plot_type}" if plot_type != "heatmap" else ""
    fig_filename_png = (
        figures_dir
        / f"{project}.{modality}.{regression_type}.celltype_similarity.{target_variable}.{effect_column}{space_suffix}{plot_suffix}.png"
    )
    fig_filename_svg = (
        figures_dir
        / f"{project}.{modality}.{regression_type}.celltype_similarity.{target_variable}.{effect_column}{space_suffix}{plot_suffix}.svg"
    )
    logger.info(f"Generating clustered {plot_type}.")

    cbar_pos = (1.05, 0.15, 0.03, 0.3) if plot_type == "dotplot" else (1.05, 0.2, 0.03, 0.6)

    sns.set_theme(style="white")
    g = sns.clustermap(
        corr_matrix,
        cmap="vlag",
        annot=True if plot_type == "heatmap" else False,
        annot_kws={"size": 8},
        fmt=".2f",
        figsize=(10, 10),
        cbar_pos=cbar_pos,
        cbar_kws={"label": "Spearman Correlation"} if plot_type == "dotplot" else None,
        vmin=-1,
        vmax=1,
        dendrogram_ratio=0.01,
    )
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    if plot_type == "dotplot":
        row_order = g.dendrogram_row.reordered_ind
        col_order = g.dendrogram_col.reordered_ind

        reordered_corr = corr_matrix.iloc[row_order, col_order]
        reordered_counts = count_matrix.iloc[row_order, col_order]

        g.ax_heatmap.clear()

        nrows, ncols = reordered_corr.shape
        x, y = np.meshgrid(np.arange(ncols), np.arange(nrows))
        x_flat = x.flatten() + 0.5
        y_flat = y.flatten() + 0.5
        c_flat = reordered_corr.values.flatten()
        counts_flat = reordered_counts.values.flatten()

        c_max = counts_flat.max() if len(counts_flat) > 0 else 0
        d_max = 20.0
        d_min = 4.0

        diameters = np.zeros_like(counts_flat, dtype=float)
        if c_max > 0:
            nonzero = counts_flat > 0
            diameters[nonzero] = np.maximum(d_min, d_max * (counts_flat[nonzero] / c_max))
        s_flat = diameters ** 2

        # Draw subtle grid lines behind dots
        g.ax_heatmap.set_xticks(np.arange(ncols + 1), minor=True)
        g.ax_heatmap.set_yticks(np.arange(nrows + 1), minor=True)
        g.ax_heatmap.grid(which="minor", color="#e0e0e0", linestyle="-", linewidth=0.5)

        mask = counts_flat > 0
        if np.any(mask):
            g.ax_heatmap.scatter(
                x_flat[mask],
                y_flat[mask],
                c=c_flat[mask],
                s=s_flat[mask],
                cmap="vlag",
                vmin=-1,
                vmax=1,
                edgecolors="#555555",
                linewidth=0.5,
            )

        g.ax_heatmap.set_xticks(np.arange(ncols) + 0.5)
        g.ax_heatmap.set_yticks(np.arange(nrows) + 0.5)
        g.ax_heatmap.set_xticklabels(reordered_corr.columns, rotation=90)
        g.ax_heatmap.set_yticklabels(reordered_corr.index, rotation=0)
        g.ax_heatmap.set_xlim(0, ncols)
        g.ax_heatmap.set_ylim(nrows, 0)

        # Legend for dot sizes anchored above colorbar
        if c_max > 0:
            nonzero_counts = counts_flat[counts_flat > 0]
            if len(nonzero_counts) > 0:
                c_min = nonzero_counts.min()
                legend_vals = np.unique(np.linspace(c_min, c_max, num=4, dtype=int))
                legend_handles = []
                for val in legend_vals:
                    d = max(d_min, d_max * (val / c_max))
                    legend_handles.append(
                        plt.scatter([], [], s=d**2, c="grey", edgecolors="#555555", linewidth=0.5)
                    )
                leg = g.ax_cbar.legend(
                    legend_handles,
                    [str(v) for v in legend_vals],
                    title="Shared Features",
                    bbox_to_anchor=(0.5, 1.15),
                    loc="lower center",
                    frameon=False,
                    labelspacing=1.8,
                    handletextpad=1.2,
                    borderpad=0.5,
                    scatterpoints=1,
                )
                leg.get_title().set_fontsize(9)
                leg.get_title().set_fontweight("bold")
                for t in leg.get_texts():
                    t.set_fontsize(8)

    title = f"Cell-Type Similarity\nModality: {modality.upper()}, Target: {target_variable.upper()}, Effect: {effect_column}"
    if feature_space != "global":
        title += f", Space: {feature_space.capitalize()}"
    if plot_type != "heatmap":
        title += f", Type: {plot_type.capitalize()}"
    g.ax_heatmap.set_title(title, pad=20)
    g.ax_heatmap.set_xlabel("Cell Type")
    g.ax_heatmap.set_ylabel("Cell Type")

    g.savefig(fig_filename_png, dpi=300, bbox_inches="tight")
    g.savefig(fig_filename_svg, dpi=300, bbox_inches="tight")
    plt.close()

    logger.info(f"Saved figures to {fig_filename_png} and {fig_filename_svg}")


if __name__ == "__main__":
    main()
