import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import seaborn as sns
import scanpy as sc
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import fdrcorrection
from cnmf import cNMF

from pseudobulk_convert import MODAL_TYPES_DICT

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare age-associated latent factors across cell types and modalities."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument(
        "--density-threshold",
        type=float,
        default=0.1,
        help="Density threshold used in cNMF.",
    )
    parser.add_argument(
        "--per-cell",
        action="store_true",
        help="Compare factor usage at the single-cell level. Default is to pseudobulk by sample_id first.",
    )
    parser.add_argument("--debug", action="store_true")
    parser.add_argument(
        "--cnmf-dir-name",
        type=str,
        default="cnmf",
        help="Name of the cnmf output directory within the latents path.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for inclusion (nominally significant factors).",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for heatmap annotation.",
    )
    parser.add_argument(
        "--corr-method",
        type=str,
        choices=["spearman", "pearson"],
        default="spearman",
        help="Correlation method to use.",
    )
    return parser.parse_args()


def load_sig_factors(results_dir, project, modality, p_threshold, per_cell):
    regression_type = "lmm" if per_cell else "pb_wls"
    fdr_file = (
        results_dir
        / "latents"
        / "cnmf"
        / ("lmm_results" if per_cell else "wls_pb_results")
        / f"{project}_combined_{modality}_{regression_type}_fdr.csv"
    )

    if not fdr_file.exists():
        logger.warning(f"File not found: {fdr_file}")
        return pd.DataFrame()

    df = pd.read_csv(fdr_file)
    if "pval_age" not in df.columns:
        logger.warning(f"'pval_age' column missing in {fdr_file}")
        return pd.DataFrame()

    sig_df = df[df["pval_age"] <= p_threshold].copy()
    sig_df["modality"] = modality
    logger.info(
        f"Loaded {len(sig_df)} nominally significant factors (pval_age <= {p_threshold}) for {modality}."
    )
    return sig_df


def main():
    args = parse_args()
    debug = args.debug

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"
    figs_dir = work_dir / "figures"

    logs_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    file_suffix = "lmm" if args.per_cell else "pb_wls"
    cnmf_dir = results_dir / "latents" / args.cnmf_dir_name

    log_filename = logs_dir / f"{args.project}_compare_latent_factors_{file_suffix}.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
    )

    logger.info(f"Command line: {' '.join(sys.argv)}")

    # 1. Identify nominally significant factors
    rna_sig = load_sig_factors(
        results_dir, args.project, "rna", args.p_threshold, args.per_cell
    )
    atac_sig = load_sig_factors(
        results_dir, args.project, "atac", args.p_threshold, args.per_cell
    )

    sig_factors = pd.concat([rna_sig, atac_sig], ignore_index=True)
    if sig_factors.empty:
        logger.error("No significant factors found.")
        sys.exit(1)

    logger.info(f"Total significant factors to compare: {len(sig_factors)}")

    annot_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"
    logger.info(f"Loading annotated data obs from {annot_file}...")
    try:
        adata = sc.read_h5ad(annot_file, backed="r")
        # Only load the columns we need to save memory
        obs_df = adata.obs[["modality", "cell_label", "sample_id"]].copy()
        # Create a safe cell label to match cNMF directory names and results
        obs_df["safe_cell_label"] = (
            obs_df["cell_label"].str.replace(" ", "_").str.replace("/", "-")
        )
    except Exception as e:
        logger.error(f"Failed to load annotated anndata object: {e}")
        sys.exit(1)

    # Dictionary to store usage dataframes for each significant factor
    # Key: {cell_type}_{modality}_K{k}_F{factor}
    factor_usages = {}
    factor_fdr_age = {}
    factor_spectra = {"rna": {}, "atac": {}}

    # Group by cell_type and modality to minimize loading cNMF objects
    grouped_factors = sig_factors.groupby(["cell_type", "modality", "k"])

    for (ct, mod, k), group in grouped_factors:
        run_name = f"{args.project}_{ct}_{mod}"

        logger.info(f"Loading cNMF usages for {run_name} K={k}")

        cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=run_name)
        try:
            cnmf_res = cnmf_obj.load_results(
                K=k, density_threshold=args.density_threshold
            )
            usage = cnmf_res[0]
            spectra = cnmf_res[1]
        except Exception as e:
            logger.error(
                f"Failed to load cNMF results for {run_name} K={k}. Error: {e}"
            )
            continue

        modality_list = MODAL_TYPES_DICT.get(mod, [mod])
        ct_obs = obs_df[
            (obs_df.modality.isin(modality_list)) & (obs_df["safe_cell_label"] == ct)
        ].copy()

        # Ensure sample_id is string
        ct_obs["sample_id"] = ct_obs["sample_id"].astype(str)

        # Filter usage index to only cells in ct_obs
        usage = usage.reindex(ct_obs.index).dropna()
        if len(usage) == 0:
            logger.warning(
                f"No cells left after intersecting cNMF usages with obs for {run_name}"
            )
            continue

        # Add sample_id for pseudobulking
        usage_with_sample = usage.merge(
            ct_obs[["sample_id"]], left_index=True, right_index=True, how="left"
        )

        if not args.per_cell:
            # Pseudobulk by taking the mean usage per sample
            # Exclude cells that might not have a sample_id
            usage_with_sample = usage_with_sample.dropna(subset=["sample_id"])
            agg_usage = usage_with_sample.groupby("sample_id").mean()
        else:
            # Per-cell comparison (warning: very memory intensive if combining across many cell types)
            # Maybe use random sampling or just use it as is?
            # Actually, per-cell comparison across DIFFERENT cell types doesn't make sense since they have different cells.
            # So if per-cell, we can only correlate within the same cell type or modality?
            # The prompt says "compare any latent factor... in other cell-types and acrouss modalities".
            # Cross-cell-type correlation REQUIRES pseudobulking to a common set of sample IDs!
            logger.warning(
                "Comparing across cell types requires a common index (like sample_id). Forcing pseudobulk for cross-comparison."
            )
            usage_with_sample = usage_with_sample.dropna(subset=["sample_id"])
            agg_usage = usage_with_sample.groupby("sample_id").mean()

        for factor_idx in group["factor"].values:
            factor_str = str(factor_idx)
            factor_int = int(factor_idx)
            fdr_age_val = group[group["factor"] == factor_idx].iloc[0].get("fdr_age", 1.0)
            
            # Extract usage and spectra based on string index match
            if factor_str in agg_usage.columns:
                col_name = f"{ct}|{mod.upper()}|K{k}|F{factor_str}"
                factor_usages[col_name] = agg_usage[factor_str]
                factor_fdr_age[col_name] = fdr_age_val
                
                # Spectra handling
                if factor_str in spectra.columns:
                    factor_spectra[mod][col_name] = spectra[factor_str]
                elif factor_int in spectra.columns:
                    factor_spectra[mod][col_name] = spectra[factor_int]
            
            # Extract usage and spectra based on integer index match
            elif factor_int in agg_usage.columns:
                col_name = f"{ct}|{mod.upper()}|K{k}|F{factor_idx}"
                factor_usages[col_name] = agg_usage[factor_int]
                factor_fdr_age[col_name] = fdr_age_val
                
                # Spectra handling
                if factor_int in spectra.columns:
                    factor_spectra[mod][col_name] = spectra[factor_int]
                elif factor_str in spectra.columns:
                    factor_spectra[mod][col_name] = spectra[factor_str]
            else:
                logger.warning(
                    f"Factor {factor_idx} not found in usage columns for {run_name} K={k}"
                )

    if not factor_usages:
        logger.error("No usage data could be extracted for any significant factors.")
        sys.exit(1)

    # Combine into a single dataframe where rows are sample_ids and columns are factors
    combined_usages = pd.DataFrame(factor_usages)

    logger.info(
        f"Combined usage matrix shape: {combined_usages.shape} (Samples x Factors)"
    )

    # Drop samples that have NaN across ALL factors?
    # Or just use pairwise complete observations for correlation.
    corr_matrix = combined_usages.corr(method=args.corr_method)

    # Fill NA correlations (e.g., if a factor has 0 variance) with 0
    corr_matrix = corr_matrix.fillna(0)

    # Calculate p-values
    logger.info("Calculating pairwise p-values and FDR...")
    cols = combined_usages.columns
    n_cols = len(cols)
    pval_matrix = pd.DataFrame(np.ones((n_cols, n_cols)), index=cols, columns=cols)
    corr_func = spearmanr if args.corr_method == "spearman" else pearsonr

    for i in range(n_cols):
        for j in range(i + 1, n_cols):
            valid_mask = (
                combined_usages.iloc[:, i].notna() & combined_usages.iloc[:, j].notna()
            )
            mask_arr = valid_mask.values
            if mask_arr.sum() >= 3:
                _, p = corr_func(
                    combined_usages.iloc[mask_arr, i], combined_usages.iloc[mask_arr, j]
                )
                p = 1.0 if np.isnan(p) else p
                pval_matrix.iloc[i, j] = p
                pval_matrix.iloc[j, i] = p

    # Extract upper triangle p-values
    upper_tri_indices = np.triu_indices(n_cols, k=1)
    pvals_upper = pval_matrix.values[upper_tri_indices]

    # Calculate FDR
    _, fdr_upper = fdrcorrection(pvals_upper)

    # Reconstruct FDR matrix
    fdr_matrix = pd.DataFrame(np.ones((n_cols, n_cols)), index=cols, columns=cols)
    fdr_matrix.values[upper_tri_indices] = fdr_upper
    fdr_matrix.values.T[upper_tri_indices] = fdr_matrix.values[upper_tri_indices]

    # Create annotation matrix for heatmap
    annot_matrix = pd.DataFrame("", index=cols, columns=cols)
    annot_matrix[fdr_matrix <= args.fdr_threshold] = "*"
    np.fill_diagonal(annot_matrix.values, "")

    # Create pairwise summary table
    pairwise_results = []
    for i in range(len(upper_tri_indices[0])):
        row_idx = upper_tri_indices[0][i]
        col_idx = upper_tri_indices[1][i]
        pairwise_results.append({
            "factor1": cols[row_idx],
            "factor2": cols[col_idx],
            "correlation": corr_matrix.iloc[row_idx, col_idx],
            "p_value": pval_matrix.iloc[row_idx, col_idx],
            "fdr": fdr_upper[i]
        })
    pairwise_df = pd.DataFrame(pairwise_results)

    # Save correlation matrix
    corr_out = (
        results_dir
        / "latents"
        / "cnmf"
        / f"{args.project}_age_factors_correlation_{file_suffix}.csv"
    )
    corr_matrix.to_csv(corr_out)

    # Save FDR matrix
    fdr_out = (
        results_dir
        / "latents"
        / "cnmf"
        / f"{args.project}_age_factors_fdr_{file_suffix}.csv"
    )
    fdr_matrix.to_csv(fdr_out)

    # Save pairwise FDR table
    fdr_table_out = (
        results_dir
        / "latents"
        / "cnmf"
        / f"{args.project}_age_factors_pairwise_fdr_table_{file_suffix}.csv"
    )
    pairwise_df.to_csv(fdr_table_out, index=False)

    logger.info(f"Saved correlation, FDR matrices, and pairwise table to {corr_out.parent}")

    # Generate clustered heatmap
    logger.info("Generating clustered heatmap...")

    # Optional: create row/col colors based on modality
    modality_colors = {"RNA": "purple", "ATAC": "green"}
    col_colors = [modality_colors[col.split("|")[1]] for col in corr_matrix.columns]

    try:
        # Increase figure size depending on number of factors
        n_factors = corr_matrix.shape[0]
        fig_size = max(10, n_factors * 0.3)

        g = sns.clustermap(
            corr_matrix,
            cmap="RdBu_r",
            vmin=-1,
            vmax=1,
            center=0,
            figsize=(fig_size, fig_size),
            xticklabels=True,
            yticklabels=True,
            col_colors=col_colors,
            row_colors=col_colors,
            linewidths=0.5,
            cbar_kws={"label": f"{args.corr_method.capitalize()} Correlation"},
            annot=annot_matrix,
            fmt="",
            annot_kws={"color": "black"},
        )

        g.ax_heatmap.set_xticklabels(
            g.ax_heatmap.get_xmajorticklabels(), fontsize=8, rotation=90
        )
        g.ax_heatmap.set_yticklabels(
            g.ax_heatmap.get_ymajorticklabels(), fontsize=8, rotation=0
        )

        # Annotate significant factors on the color bands
        for tick in g.ax_heatmap.get_xticklabels():
            col_name = tick.get_text()
            if factor_fdr_age.get(col_name, 1.0) <= args.fdr_threshold:
                x = tick.get_position()[0]
                g.ax_col_colors.text(x, 0.5, '+', ha='center', va='center', color='white', fontweight='bold', fontsize=8)

        for tick in g.ax_heatmap.get_yticklabels():
            row_name = tick.get_text()
            if factor_fdr_age.get(row_name, 1.0) <= args.fdr_threshold:
                y = tick.get_position()[1]
                g.ax_row_colors.text(0.5, y, '+', ha='center', va='center', color='white', fontweight='bold', fontsize=8, rotation=90)

        heatmap_out_png = (
            figs_dir / f"{args.project}_age_factors_clustermap_{file_suffix}.png"
        )
        heatmap_out_svg = (
            figs_dir / f"{args.project}_age_factors_clustermap_{file_suffix}.svg"
        )

        g.savefig(heatmap_out_png, dpi=300, bbox_inches="tight")
        g.savefig(heatmap_out_svg, bbox_inches="tight")
        logger.info(f"Saved clustered heatmap to {heatmap_out_png}")

        # Generate narrower heatmap containing only factors with at least 1 significant correlation
        # Exclude the diagonal which is 1.0; off-diagonal elements in fdr_matrix represent the actual comparisons
        sig_mask = (fdr_matrix <= args.fdr_threshold).any(axis=0)
        sig_factors_cols = fdr_matrix.columns[sig_mask]

        if len(sig_factors_cols) > 1:
            logger.info(f"Generating narrower clustered heatmap with {len(sig_factors_cols)} factors...")
            narrow_corr = corr_matrix.loc[sig_factors_cols, sig_factors_cols]
            narrow_annot = annot_matrix.loc[sig_factors_cols, sig_factors_cols]
            narrow_colors = [modality_colors[col.split("|")[1]] for col in sig_factors_cols]

            n_factors_narrow = narrow_corr.shape[0]
            fig_size_narrow = max(8, n_factors_narrow * 0.3)

            g_narrow = sns.clustermap(
                narrow_corr,
                cmap="RdBu_r",
                vmin=-1,
                vmax=1,
                center=0,
                figsize=(fig_size_narrow, fig_size_narrow),
                xticklabels=True,
                yticklabels=True,
                col_colors=narrow_colors,
                row_colors=narrow_colors,
                linewidths=0.5,
                cbar_kws={"label": f"{args.corr_method.capitalize()} Correlation"},
                annot=narrow_annot,
                fmt="",
                annot_kws={"color": "black"},
            )

            g_narrow.ax_heatmap.set_xticklabels(
                g_narrow.ax_heatmap.get_xmajorticklabels(), fontsize=8, rotation=90
            )
            g_narrow.ax_heatmap.set_yticklabels(
                g_narrow.ax_heatmap.get_ymajorticklabels(), fontsize=8, rotation=0
            )

            for tick in g_narrow.ax_heatmap.get_xticklabels():
                col_name = tick.get_text()
                if factor_fdr_age.get(col_name, 1.0) <= args.fdr_threshold:
                    x = tick.get_position()[0]
                    g_narrow.ax_col_colors.text(x, 0.5, '+', ha='center', va='center', color='white', fontweight='bold', fontsize=8)

            for tick in g_narrow.ax_heatmap.get_yticklabels():
                row_name = tick.get_text()
                if factor_fdr_age.get(row_name, 1.0) <= args.fdr_threshold:
                    y = tick.get_position()[1]
                    g_narrow.ax_row_colors.text(0.5, y, '+', ha='center', va='center', color='white', fontweight='bold', fontsize=8, rotation=90)

            narrow_heatmap_out_png = figs_dir / f"{args.project}_age_factors_clustermap_narrow_{file_suffix}.png"
            narrow_heatmap_out_svg = figs_dir / f"{args.project}_age_factors_clustermap_narrow_{file_suffix}.svg"
            g_narrow.savefig(narrow_heatmap_out_png, dpi=300, bbox_inches="tight")
            g_narrow.savefig(narrow_heatmap_out_svg, bbox_inches="tight")
            logger.info(f"Saved narrower clustered heatmap to {narrow_heatmap_out_png}")
        else:
            logger.info("Not enough factors with significant correlations to generate a narrower heatmap.")

    except Exception as e:
        logger.error(f"Failed to generate clustermap: {e}")

    # Generate within-modality spectra correlation heatmaps
    logger.info("Generating within-modality spectra correlation heatmaps...")
    for mod in ["rna", "atac"]:
        if not factor_spectra[mod]:
            logger.info(f"No spectra extracted for {mod.upper()}, skipping spectra correlation.")
            continue
            
        logger.info(f"Processing {mod.upper()} spectra correlation...")
        spectra_df = pd.DataFrame(factor_spectra[mod]).fillna(0)
        spectra_corr = spectra_df.corr(method=args.corr_method)
        
        spectra_corr_out = results_dir / "latents" / "cnmf" / f"{args.project}_age_factors_spectra_correlation_{mod}_{file_suffix}.csv"
        spectra_corr.to_csv(spectra_corr_out)
        
        mod_color = modality_colors[mod.upper()]
        spec_colors = [mod_color] * len(spectra_corr.columns)
        
        try:
            n_factors_mod = spectra_corr.shape[0]
            fig_size_mod = max(8, n_factors_mod * 0.3)
            
            g_spec = sns.clustermap(
                spectra_corr, 
                cmap="RdBu_r", 
                vmin=-1, 
                vmax=1,
                center=0,
                figsize=(fig_size_mod, fig_size_mod),
                xticklabels=True,
                yticklabels=True,
                col_colors=spec_colors,
                row_colors=spec_colors,
                linewidths=0.5,
                cbar_kws={"label": f"{args.corr_method.capitalize()} Correlation"}
            )
            
            g_spec.ax_heatmap.set_xticklabels(g_spec.ax_heatmap.get_xmajorticklabels(), fontsize=8, rotation=90)
            g_spec.ax_heatmap.set_yticklabels(g_spec.ax_heatmap.get_ymajorticklabels(), fontsize=8, rotation=0)
            
            for tick in g_spec.ax_heatmap.get_xticklabels():
                col_name = tick.get_text()
                if factor_fdr_age.get(col_name, 1.0) <= args.fdr_threshold:
                    x = tick.get_position()[0]
                    g_spec.ax_col_colors.text(x, 0.5, '+', ha='center', va='center', color='white', fontweight='bold', fontsize=8)

            for tick in g_spec.ax_heatmap.get_yticklabels():
                row_name = tick.get_text()
                if factor_fdr_age.get(row_name, 1.0) <= args.fdr_threshold:
                    y = tick.get_position()[1]
                    g_spec.ax_row_colors.text(0.5, y, '+', ha='center', va='center', color='white', fontweight='bold', fontsize=8, rotation=90)
                    
            heatmap_out_png = figs_dir / f"{args.project}_age_factors_spectra_clustermap_{mod}_{file_suffix}.png"
            heatmap_out_svg = figs_dir / f"{args.project}_age_factors_spectra_clustermap_{mod}_{file_suffix}.svg"
            g_spec.savefig(heatmap_out_png, dpi=300, bbox_inches="tight")
            g_spec.savefig(heatmap_out_svg, bbox_inches="tight")
            logger.info(f"Saved {mod.upper()} spectra clustermap to {heatmap_out_png}")
        except Exception as e:
            logger.error(f"Failed to generate spectra clustermap for {mod}: {e}")

if __name__ == "__main__":
    main()
