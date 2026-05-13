import argparse
import logging
from pathlib import Path
import sys

import pandas as pd
from statsmodels.stats import multitest as smm

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
    sig_results_df = results_df.loc[results_df["fdr_age"] <= 0.05]
    sig_results_df.to_csv(sig_out_file, index=False)
    logger.info(f"Saved significant combined FDR results to {sig_out_file}")


if __name__ == "__main__":
    main()
