import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import scanpy as sc
import statsmodels.formula.api as smf
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
        description="Run a mixed effects linear regression model between age and cNMF latent factors."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--modality", type=str, default="rna", choices=["rna", "atac"])
    parser.add_argument(
        "--cell-type", type=str, required=True, help="Target cell type."
    )
    parser.add_argument("--k", type=int, required=True, help="Selected K for cNMF.")
    parser.add_argument(
        "--density-threshold",
        type=float,
        default=0.1,
        help="Density threshold used in cNMF.",
    )
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    debug = args.debug

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    run_name = f"{args.project}_{safe_ct}_{args.modality}"

    log_filename = (
        logs_dir / f"{args.project}_{args.modality}_{safe_ct}_lmm_k{args.k}.log"
    )
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
    )

    logger.info(f"Command line: {' '.join(sys.argv)}")

    annot_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"

    # Load Annotated Data
    logger.info(f"Loading annotated data from {annot_file}...")
    try:
        # Load in backed mode since we only need the .obs metadata, saving memory and time
        adata = sc.read_h5ad(annot_file, backed="r")
    except Exception as e:
        logger.error(f"Failed to load annotated anndata object: {e}")
        sys.exit(1)

    logger.info(
        f"Subsetting metadata by modality: {args.modality} and cell type: {args.cell_type}"
    )
    modality_list = MODAL_TYPES_DICT.get(args.modality)

    try:
        # We only need the observation metadata for the mixed model
        obs_df = adata.obs[
            (adata.obs.modality.isin(modality_list))
            & (adata.obs["cell_label"] == args.cell_type)
        ]
    except Exception as e:
        logger.error(f"Failed to filter anndata observations: {e}")
        sys.exit(1)

    if len(obs_df) == 0:
        logger.error(
            f"No cells found for cell type {args.cell_type} in modality {args.modality}."
        )
        sys.exit(1)

    logger.info(f"Found {len(obs_df)} cells for {args.cell_type}.")

    cnmf_dir = results_dir / "cnmf"
    logger.info(
        f"Loading cNMF results for {run_name} with K={args.k} and density_threshold={args.density_threshold}"
    )

    cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=run_name)

    try:
        # load_results returns usage, spectra_scores, spectra_tpm, top_genes
        usage, *_ = cnmf_obj.load_results(
            K=args.k, density_threshold=args.density_threshold
        )
    except Exception as e:
        logger.error(
            f"Failed to load cNMF results. Ensure cNMF consensus has finished for K={args.k}. Error: {e}"
        )
        sys.exit(1)

    # cNMF usage index corresponds to the cell IDs from the anndata.
    metadata = obs_df[["age", "sample_id"]].copy()
    metadata["age"] = pd.to_numeric(metadata["age"], errors="coerce")

    # Drop rows missing age or sample_id
    metadata = metadata.dropna(subset=["age", "sample_id"])

    # Merge metadata with usage by cell IDs
    df = pd.merge(usage, metadata, left_index=True, right_index=True, how="inner")

    if len(df) == 0:
        logger.error(
            "No overlapping cells between cNMF results and metadata after dropping NA age/sample_id."
        )
        sys.exit(1)

    logger.info(
        f"Running Linear Mixed Effects Model (MixedLM) for {len(df)} cells across {usage.shape[1]} latent factors."
    )
    logger.info("Formula: Factor ~ age + (1|sample_id)")

    lmm_results = []
    # Usage columns represent latent factors (typically numbers or string numbers like '1', '2')
    factor_cols = usage.columns.tolist()

    for factor in factor_cols:
        logger.info(f"Fitting MixedLM for latent factor '{factor}'...")

        # Prepare dataframe for this factor (rename column to avoid formula issues with numeric column names)
        # Using a temporary column name 'latent_factor'
        formula_df = df[[factor, "age", "sample_id"]].copy()
        formula_df.rename(columns={factor: "latent_factor"}, inplace=True)

        try:
            # Fit the model
            md = smf.mixedlm(
                "latent_factor ~ age", formula_df, groups=formula_df["sample_id"]
            )
            mdf = md.fit()

            # Extract statistics for 'age'
            age_coef = mdf.params.get("age", np.nan)
            age_pval = mdf.pvalues.get("age", np.nan)
            age_se = mdf.bse.get("age", np.nan)
            converged = mdf.converged

            lmm_results.append(
                {
                    "factor": factor,
                    "coef_age": age_coef,
                    "se_age": age_se,
                    "pval_age": age_pval,
                    "converged": converged,
                }
            )

            logger.info(
                f"Factor {factor}: coef_age={age_coef:.4f}, pval_age={age_pval:.4e}, converged={converged}"
            )

        except Exception as e:
            logger.error(f"MixedLM failed for factor '{factor}': {e}")
            lmm_results.append(
                {
                    "factor": factor,
                    "coef_age": np.nan,
                    "se_age": np.nan,
                    "pval_age": np.nan,
                    "converged": False,
                }
            )

    results_df = pd.DataFrame(lmm_results)

    # Save results
    out_dir = results_dir / "lmm_results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_file = out_dir / f"{args.project}_{safe_ct}_{args.modality}_k{args.k}_lmm.csv"
    results_df.to_csv(out_file, index=False)

    logger.info(
        f"Linear Mixed Effects Model analysis complete. Results saved to:\n  {out_file}"
    )


if __name__ == "__main__":
    main()
