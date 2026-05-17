import argparse
import logging
from pathlib import Path
import sys

import pandas as pd

# Import functions from pseudobulk_regression
from pseudobulk_regression import compute_fdr

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run post-processing for mediation analysis."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--endo-modality", type=str, default="rna")
    parser.add_argument("--exog-modality", type=str, default="atac")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"

    log_filename = f"{logs_dir}/post_mediation_{args.endo_modality}_{args.exog_modality}_wls.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    # Load cis fdr filtered results to get the list of cell types
    # This ensures we only look for cell types that had significant cis pairs
    input_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.wls.cis.fdr_filtered.csv"
    )
    
    if not input_file.exists():
        logger.error(f"Significant cis correlations input file not found: {input_file}")
        sys.exit(1)
        
    sig_results = pd.read_csv(input_file)
    cell_types = sig_results["tissue"].unique()
    
    logger.info(f"Looking for uncorrected mediation results for {len(cell_types)} cell types...")

    all_results = []
    
    for cell_type in cell_types:
        ct_file = (
            results_dir
            / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{cell_type}.wls.mediation.csv"
        )
        if ct_file.exists():
            ct_df = pd.read_csv(ct_file)
            all_results.append(ct_df)
            logger.info(f"Loaded {len(ct_df)} results for {cell_type}")
        else:
            logger.warning(f"File not found for {cell_type}: {ct_file}")
            
    if not all_results:
        logger.error("No mediation results found to process.")
        sys.exit(1)

    results_df = pd.concat(all_results, ignore_index=True)
    
    # Filter out missing results (should already be done by cis_correlation_mediation script, but just in case)
    results_df = results_df.dropna(subset=["acme_estimate"])

    if len(results_df) == 0:
        logger.info("All loaded mediation runs failed or resulted in NaNs. Exiting.")
        return

    # Calculate FDR on ACME (Average Causal Mediation Effect) p-value
    results_df["acme_fdr"] = compute_fdr(results_df["acme_pval"])

    total_sig = results_df.loc[results_df["acme_fdr"] <= 0.05].shape[0]
    logger.info(
        f"Found {total_sig} significantly mediated pairs across all cell types (ACME FDR <= 0.05)"
    )

    out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.wls.mediation.csv"
    )
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved combined mediation results to {out_file}")

    sig_out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.wls.mediation.fdr_filtered.csv"
    )
    sig_results_df = results_df.loc[results_df["acme_fdr"] <= 0.05]
    sig_results_df.to_csv(sig_out_file, index=False)
    logger.info(f"Saved significant combined mediation results to {sig_out_file}")


if __name__ == "__main__":
    main()
