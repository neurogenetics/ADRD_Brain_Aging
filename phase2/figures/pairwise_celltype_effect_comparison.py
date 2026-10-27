#!/usr/bin/env python3
"""
Visualize pairwise target-variable effect comparisons between two cell-types.
Includes features in the union of significant (FDR-filtered) results across the
two cell-types and colors points by whether each feature is significant in both
cell-types or specific to one cell-type.
"""

import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare target variable effects pairwise between two cell-types."
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
        "--cell-type-1",
        type=str,
        required=True,
        help="First cell-type name (plotted on x-axis).",
    )
    parser.add_argument(
        "--cell-type-2",
        type=str,
        required=True,
        help="Second cell-type name (plotted on y-axis).",
    )
    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
        help="Regression method to use (default: 'wls').",
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
        help="Effect column to compare (default: 'coef').",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def normalize_cell_type(cell_type: str) -> str:
    return cell_type.replace(" ", "_").replace("/", "-")


def load_full_results_for_celltype(
    results_dir: Path,
    project: str,
    modality: str,
    cell_type: str,
    regression_type: str,
    target_variable: str,
) -> pd.DataFrame:
    """
    Load full regression results for a specific cell-type.
    Checks the per-cell-type result file first, then falls back to the
    aggregated all_celltypes file.
    """
    per_ct_file = (
        results_dir
        / f"{project}.{modality}.{cell_type}.{regression_type}.{target_variable}.csv"
    )
    if per_ct_file.exists():
        logger.info("Loading full results for %s from %s", cell_type, per_ct_file)
        df = pd.read_csv(per_ct_file)
        return df

    all_ct_file = (
        results_dir
        / f"{project}.all_celltypes.{modality}.{regression_type}.{target_variable}.csv"
    )
    if all_ct_file.exists():
        logger.info(
            "Per-cell-type file not found; loading %s from aggregated %s",
            cell_type,
            all_ct_file,
        )
        df = pd.read_csv(all_ct_file)
        if "tissue" in df.columns:
            df = df.loc[df["tissue"] == cell_type].copy()
        return df

    raise FileNotFoundError(
        f"Could not find full results for cell-type '{cell_type}'. "
        f"Checked: {per_ct_file} and {all_ct_file}"
    )


def load_sig_features_for_celltype(
    results_dir: Path,
    project: str,
    modality: str,
    cell_type: str,
    regression_type: str,
    target_variable: str,
) -> set[str]:
    """
    Load significant (FDR-filtered) features for a specific cell-type.
    Checks the per-cell-type FDR-filtered file first, then falls back to the
    aggregated all_celltypes FDR-filtered file.
    """
    per_ct_fdr_file = (
        results_dir
        / f"{project}.{modality}.{cell_type}.{regression_type}_fdr_filtered.{target_variable}.csv"
    )
    if per_ct_fdr_file.exists():
        logger.info(
            "Loading FDR-filtered results for %s from %s", cell_type, per_ct_fdr_file
        )
        df = pd.read_csv(per_ct_fdr_file)
        return set(df["feature"].dropna().unique())

    all_ct_fdr_file = (
        results_dir
        / f"{project}.{modality}.all_celltypes.{regression_type}_fdr_filtered.{target_variable}.csv"
    )
    if all_ct_fdr_file.exists():
        df = pd.read_csv(all_ct_fdr_file)
        if "tissue" in df.columns and (df["tissue"] == cell_type).any():
            logger.info(
                "Loading FDR-filtered results for %s from aggregated %s",
                cell_type,
                all_ct_fdr_file,
            )
            sub_df = df.loc[df["tissue"] == cell_type]
            return set(sub_df["feature"].dropna().unique())

    raise FileNotFoundError(
        f"Could not find FDR-filtered results for cell-type '{cell_type}'. "
        f"Checked: {per_ct_fdr_file} and {all_ct_fdr_file} (tissue == '{cell_type}')"
    )


def main():
    args = parse_args()
    debug = args.debug

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figures_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"
    figures_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    ct1 = normalize_cell_type(args.cell_type_1)
    ct2 = normalize_cell_type(args.cell_type_2)
    modality = args.modality.lower()
    regression_type = args.regression_type
    target_variable = args.target_variable
    effect_column = args.effect_column
    project = args.project

    log_filename = (
        logs_dir
        / f"{modality}_{regression_type}_{ct1}_vs_{ct2}_{target_variable}_pairwise_comparison.log"
    )
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logger.info("Command line: %s", " ".join(sys.argv))
    logger.info("Logging configured. Writing to %s", log_filename)

    # 1. Load significant features for both cell-types
    try:
        sig_1 = load_sig_features_for_celltype(
            results_dir, project, modality, ct1, regression_type, target_variable
        )
        sig_2 = load_sig_features_for_celltype(
            results_dir, project, modality, ct2, regression_type, target_variable
        )
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    union_sig = sig_1.union(sig_2)
    shared_sig = sig_1.intersection(sig_2)
    logger.info("%s significant features: %d", ct1, len(sig_1))
    logger.info("%s significant features: %d", ct2, len(sig_2))
    logger.info(
        "Union of significant features: %d (Shared in both: %d)",
        len(union_sig),
        len(shared_sig),
    )

    if len(union_sig) == 0:
        logger.warning(
            "No significant features found in either %s or %s. Skipping plot.",
            ct1,
            ct2,
        )
        return

    # 2. Load full results for both cell-types
    try:
        full_1 = load_full_results_for_celltype(
            results_dir, project, modality, ct1, regression_type, target_variable
        )
        full_2 = load_full_results_for_celltype(
            results_dir, project, modality, ct2, regression_type, target_variable
        )
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    for label, df in [(ct1, full_1), (ct2, full_2)]:
        if "feature" not in df.columns or effect_column not in df.columns:
            logger.error(
                "Required columns ('feature', '%s') missing in full results for %s.",
                effect_column,
                label,
            )
            sys.exit(1)

    # Filter to union of significant features and align on feature
    sub_1 = (
        full_1.loc[full_1["feature"].isin(union_sig), ["feature", effect_column]]
        .drop_duplicates(subset=["feature"])
        .rename(columns={effect_column: f"{effect_column}_{ct1}"})
    )
    sub_2 = (
        full_2.loc[full_2["feature"].isin(union_sig), ["feature", effect_column]]
        .drop_duplicates(subset=["feature"])
        .rename(columns={effect_column: f"{effect_column}_{ct2}"})
    )

    merged = pd.merge(sub_1, sub_2, on="feature", how="inner").dropna(
        subset=[f"{effect_column}_{ct1}", f"{effect_column}_{ct2}"]
    )

    if merged.empty:
        logger.warning(
            "No overlapping features with valid '%s' values between %s and %s.",
            effect_column,
            ct1,
            ct2,
        )
        return

    # Classify each feature by significance membership
    in_1 = merged["feature"].isin(sig_1)
    in_2 = merged["feature"].isin(sig_2)

    cond_both = in_1 & in_2
    cond_1_only = in_1 & ~in_2
    cond_2_only = ~in_1 & in_2

    n_both = int(cond_both.sum())
    n_1_only = int(cond_1_only.sum())
    n_2_only = int(cond_2_only.sum())

    label_both = f"Both ({n_both})"
    label_1_only = f"{ct1} only ({n_1_only})"
    label_2_only = f"{ct2} only ({n_2_only})"

    merged["sig_group"] = np.select(
        [cond_both, cond_1_only, cond_2_only],
        [label_both, label_1_only, label_2_only],
        default="Other",
    )

    # Order plotting so "Both" is drawn on top
    plot_order = [label_1_only, label_2_only, label_both]
    merged["sig_group"] = pd.Categorical(
        merged["sig_group"], categories=plot_order, ordered=True
    )
    merged = merged.sort_values("sig_group")

    # Compute Spearman correlation across the union
    x_col = f"{effect_column}_{ct1}"
    y_col = f"{effect_column}_{ct2}"
    if len(merged) >= 2:
        spearman_r, spearman_p = stats.spearmanr(merged[x_col], merged[y_col])
    else:
        spearman_r, spearman_p = np.nan, np.nan

    logger.info(
        "Union features plotted: %d | Spearman r = %.4f (p = %.3e)",
        len(merged),
        spearman_r,
        spearman_p,
    )

    if debug:
        print(merged.head(10))

    # 3. Generate Scatter Plot
    sns.set_theme(style="whitegrid")
    fig, ax = plt.subplots(figsize=(8, 8))

    palette = {
        label_1_only: "#1f77b4",
        label_2_only: "#ff7f0e",
        label_both: "#6a3d9a",
    }

    sns.scatterplot(
        data=merged,
        x=x_col,
        y=y_col,
        hue="sig_group",
        hue_order=plot_order,
        palette=palette,
        alpha=0.75,
        s=45,
        edgecolor="none",
        ax=ax,
    )

    # Zero reference lines
    ax.axhline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8, alpha=0.6)

    p_str = f"{spearman_p:.2e}" if spearman_p > 0 else "< 1e-300"
    ax.set_title(
        f"{ct1} vs {ct2} — {target_variable.upper()} Effect ({modality.upper()}, {regression_type.upper()})\n"
        f"Union N = {len(merged)} | Spearman r = {spearman_r:.3f} (p = {p_str})",
        fontsize=12,
    )
    ax.set_xlabel(f"{ct1} Effect ({effect_column})", fontsize=11)
    ax.set_ylabel(f"{ct2} Effect ({effect_column})", fontsize=11)
    ax.legend(title="Significance (FDR < 0.05)", loc="best", frameon=True)

    plt.tight_layout()

    fig_png = (
        figures_dir
        / f"{project}.{modality}.{regression_type}.pairwise_effect.{ct1}_vs_{ct2}.{target_variable}.{effect_column}.png"
    )
    fig_svg = (
        figures_dir
        / f"{project}.{modality}.{regression_type}.pairwise_effect.{ct1}_vs_{ct2}.{target_variable}.{effect_column}.svg"
    )

    plt.savefig(fig_png, dpi=300, bbox_inches="tight")
    plt.savefig(fig_svg, bbox_inches="tight")
    plt.close()

    logger.info("Saved pairwise comparison figure (PNG) to %s", fig_png)
    logger.info("Saved pairwise comparison figure (SVG) to %s", fig_svg)


if __name__ == "__main__":
    main()
