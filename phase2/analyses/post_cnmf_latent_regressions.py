import argparse
import logging
from pathlib import Path
import sys

import pandas as pd
import numpy as np
from statsmodels.stats import multitest as smm
from kneed import KneeLocator
import matplotlib.pyplot as plt
from adjustText import adjust_text

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def compute_fdr(pvalues):
    """Compute Benjamini-Hochberg FDR."""
    # Ensure there are no NaNs in pvalues before passing to fdrcorrection
    mask = pvalues.notna()
    fdr = pd.Series(index=pvalues.index, dtype=float)
    if mask.sum() > 0:
        bh_adj = smm.fdrcorrection(pvalues[mask])
        fdr[mask] = bh_adj[1]
    return fdr


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run post-processing for cNMF latent factor linear mixed effects models."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--modality", type=str, default="rna", choices=["rna", "atac"])
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"

    log_filename = f"{logs_dir}/post_cnmf_latent_regressions_{args.modality}.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    cnmf_dir = results_dir / "latents" / "cnmf"
    lmm_results_dir = cnmf_dir / "lmm_results"

    if not lmm_results_dir.exists():
        logger.error(f"Results directory does not exist: {lmm_results_dir}")
        sys.exit(1)

    logger.info(f"Searching for {args.modality} lmm results in {lmm_results_dir}")

    # File pattern: {project}_{safe_ct}_{modality}_lmm.csv
    pattern = f"{args.project}_*_{args.modality}_lmm.csv"
    file_paths = list(lmm_results_dir.glob(pattern))

    if not file_paths:
        logger.error(f"No files matching pattern {pattern} found in {lmm_results_dir}")
        sys.exit(1)

    combined_data = []
    prefix = f"{args.project}_"
    suffix = f"_{args.modality}_lmm.csv"

    for fp in file_paths:
        try:
            filename = fp.name
            if filename.startswith(prefix) and filename.endswith(suffix):
                # Extract the cell type part
                safe_ct = filename[len(prefix) : -len(suffix)]
            else:
                logger.warning(
                    f"File {filename} does not match expected naming convention, skipping extraction."
                )
                continue

            df = pd.read_csv(fp)
            df["cell_type"] = safe_ct
            combined_data.append(df)
            logger.debug(f"Loaded {len(df)} factors for {safe_ct}")
        except Exception as e:
            logger.error(f"Error reading {fp}: {e}")

    if not combined_data:
        logger.error("No valid data could be combined.")
        sys.exit(1)

    results_df = pd.concat(combined_data, ignore_index=True)
    logger.info(
        f"Combined {len(combined_data)} cell-type results, total rows: {len(results_df)}"
    )

    # Compute FDR for 'pval_age' across all combined results
    if "pval_age" not in results_df.columns:
        logger.error("Column 'pval_age' not found in combined data. Cannot compute FDR.")
        sys.exit(1)

    # Some regressions might have failed and yielded NaN for pval_age.
    # We only compute FDR on valid p-values.
    results_df["fdr_age"] = compute_fdr(results_df["pval_age"])

    # Log nominally significant results
    nominal_sig_df = results_df.loc[results_df["pval_age"] <= 0.05]
    total_nominal_sig = nominal_sig_df.shape[0]
    logger.info(
        f"Found {total_nominal_sig} nominally significant factors across all cell types (pval_age <= 0.05)"
    )
    if total_nominal_sig > 0:
        logger.info(f"Nominally significant dataframe subset:\n{nominal_sig_df.to_string(index=False)}")

    total_sig = results_df.loc[results_df["fdr_age"] <= 0.05].shape[0]
    logger.info(
        f"Found {total_sig} significant factors across all cell types (FDR_age <= 0.05)"
    )

    out_file = lmm_results_dir / f"{args.project}_combined_{args.modality}_lmm_fdr.csv"
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved combined FDR results to {out_file}")

    sig_out_file = (
        lmm_results_dir / f"{args.project}_combined_{args.modality}_lmm_fdr_filtered.csv"
    )
    sig_results_df = results_df.loc[results_df["fdr_age"] <= 0.05].copy()
    sig_results_df.to_csv(sig_out_file, index=False)
    logger.info(f"Saved significant combined FDR results to {sig_out_file}")

    # Extract top features for significant factors
    if total_sig > 0:
        logger.info("Extracting top features using the elbow method for significant factors...")
        top_features_list = []
        
        for row in sig_results_df.itertuples():
            ct = row.cell_type
            k = row.k
            factor = str(row.factor)
            
            run_name = f"{args.project}_{ct}_{args.modality}"
            run_dir = cnmf_dir / run_name
            
            if not run_dir.exists():
                logger.warning(f"cNMF run directory not found: {run_dir}. Skipping feature extraction for {ct} K={k} Factor={factor}.")
                continue
                
            # Find the spectra score file
            score_files = list(run_dir.glob(f"*spectra_score.k_{k}.dt_*.txt"))
            if not score_files:
                logger.warning(f"No spectra score file found for K={k} in {run_dir}. Skipping.")
                continue
            
            score_file = score_files[0]
            try:
                scores_df = pd.read_csv(score_file, sep='\t', index_col=0)
            except Exception as e:
                logger.error(f"Failed to read {score_file}: {e}")
                continue
                
            # cNMF index contains factors, columns are features
            # Ensure the factor exists in the index (might be int or str in index)
            factor_idx = None
            for idx in scores_df.index:
                if str(idx) == factor:
                    factor_idx = idx
                    break
                    
            if factor_idx is None:
                logger.warning(f"Factor {factor} not found in index of {score_file}. Skipping.")
                continue
                
            # Extract scores for the factor
            factor_scores = scores_df.loc[factor_idx]
            
            # Sort scores descending
            sorted_scores = factor_scores.sort_values(ascending=False)
            
            # Filter to positive scores (optional but recommended for cNMF)
            positive_scores = sorted_scores[sorted_scores > 0]
            
            if len(positive_scores) < 3:
                logger.warning(f"Not enough positive scores for {ct} K={k} Factor={factor}. Taking top 10 as fallback.")
                top_features = sorted_scores.head(10)
            else:
                y = positive_scores.values
                x = np.arange(len(y))
                
                try:
                    kneedle = KneeLocator(x, y, S=1.0, curve='convex', direction='decreasing')
                    elbow_idx = kneedle.knee
                    
                    if elbow_idx is not None:
                        # Extract features up to the elbow point (inclusive)
                        top_features = positive_scores.iloc[:elbow_idx + 1]
                    else:
                        logger.warning(f"Elbow not found for {ct} K={k} Factor={factor}. Fallback to top 100.")
                        top_features = sorted_scores.head(100)
                except Exception as e:
                    logger.warning(f"KneeLocator failed for {ct} K={k} Factor={factor}: {e}. Fallback to top 100.")
                    top_features = sorted_scores.head(100)
                    
            for feat, score in top_features.items():
                top_features_list.append({
                    "cell_type": ct,
                    "k": k,
                    "factor": factor,
                    "feature": feat,
                    "score": score
                })
                
            # Create a visualization of the features
            try:
                fig, ax = plt.subplots(figsize=(10, 6))
                
                # Get all feature scores, sorted
                all_y = sorted_scores.values
                all_x = np.arange(len(all_y))
                all_features = sorted_scores.index
                
                num_top = len(top_features)
                
                # Plot non-selected features in grey
                ax.scatter(all_x[num_top:], all_y[num_top:], color='lightgray', s=10, alpha=0.6, label='Other Features')
                
                # Plot selected top features in red
                ax.scatter(all_x[:num_top], all_y[:num_top], color='red', s=15, alpha=0.8, label='Selected Top Features')
                
                # Label the top 5 features using adjustText to prevent overlap
                texts = []
                for i in range(min(5, len(all_y))):
                    texts.append(
                        ax.text(
                            all_x[i], 
                            all_y[i], 
                            all_features[i], 
                            fontsize=9, 
                            color='black'
                        )
                    )
                
                if texts:
                    try:
                        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
                    except Exception as e:
                        logger.warning(f"adjust_text failed: {e}")
                
                # Optionally add a line at the elbow point
                if num_top > 0 and num_top < len(all_y):
                    ax.axvline(x=num_top-1, color='black', linestyle='--', alpha=0.5, label='Elbow cutoff')
                
                ax.set_title(f"cNMF Feature Scores\nCell Type: {ct} | Modality: {args.modality.upper()} | K: {k} | Factor: {factor}")
                ax.set_xlabel("Feature Rank")
                ax.set_ylabel("Spectra Score")
                ax.legend()
                
                # Save figures
                figs_dir = work_dir / "figures"
                figs_dir.mkdir(parents=True, exist_ok=True)
                base_fig_path = figs_dir / f"{args.project}_{ct}_{args.modality}_k{k}_factor{factor}_scores"
                
                fig.savefig(f"{base_fig_path}.png", bbox_inches='tight', dpi=300)
                fig.savefig(f"{base_fig_path}.svg", bbox_inches='tight')
                plt.close(fig)
                logger.debug(f"Saved visualization to {base_fig_path}.png/.svg")
                
            except Exception as e:
                logger.error(f"Failed to generate visualization for {ct} K={k} Factor={factor}: {e}")
                
        if top_features_list:
            top_features_df = pd.DataFrame(top_features_list)
            top_features_out = lmm_results_dir / f"{args.project}_combined_{args.modality}_lmm_top_features.csv"
            top_features_df.to_csv(top_features_out, index=False)
            logger.info(f"Saved extracted top features to {top_features_out}")
        else:
            logger.info("No top features were successfully extracted.")


if __name__ == "__main__":
    main()
