#!/usr/bin/env python3
"""
Model the correlation between cell counts and disease target variable across all cell types to check for WLS weight bias.
"""

import sys
import logging
import argparse
from pathlib import Path
import pandas as pd
import statsmodels.formula.api as smf

# Configure logging
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Model the correlation between cell counts and target variable across all cell types to check for WLS weight bias."
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
        "--target-variable",
        type=str,
        default="dx",
        help="The primary disease target variable column in covariates (default: 'dx').",
    )
    parser.add_argument(
        "--weight-term",
        type=str,
        default="cell_counts",
        help="The weight term being evaluated (e.g. cell_counts)",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    info_dir = work_dir / "sample_info"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"

    # Configure logging to file and stdout
    log_filename = logs_dir / f"{args.project}_batch_cell_counts_regression_disease.log"
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logger.info("Command line: %s", " ".join(sys.argv))
    logger.info("Logging configured. Writing to %s", log_filename)

    results_dir.mkdir(parents=True, exist_ok=True)

    project = args.project
    weight_term = args.weight_term
    target_variable = args.target_variable

    # Find all final_covariates files for the project
    search_pattern = f"{project}.*.final_covariates.csv"
    covariates_files = list(info_dir.glob(search_pattern))

    if not covariates_files:
        logger.error(
            "No final_covariates.csv files found in %s matching %s",
            info_dir,
            search_pattern,
        )
        sys.exit(1)

    logger.info("Found %d covariate files to process.", len(covariates_files))

    all_results = []

    for cov_file in covariates_files:
        # Filename format: project.cell_type.modality.final_covariates.csv
        parts = cov_file.name.split(".")
        if len(parts) < 5:
            logger.warning(
                "Skipping %s: does not match expected naming convention.", cov_file.name
            )
            continue

        modality = parts[-3]
        cell_type = ".".join(parts[1:-3])

        logger.info("Processing cell_type: %s, modality: %s", cell_type, modality)

        covars_df = pd.read_csv(cov_file, index_col=0)

        if weight_term not in covars_df.columns:
            logger.warning(
                "Weight term '%s' not found in %s. Skipping.",
                weight_term,
                cov_file.name,
            )
            continue

        if target_variable not in covars_df.columns:
            logger.warning(
                "Target variable '%s' not found in %s. Skipping.",
                target_variable,
                cov_file.name,
            )
            continue

        # Identify PCA terms, limiting to the first 4 to match regression modeling
        pca_terms = [col for col in covars_df.columns if col.startswith("PCA_")]
        pca_terms = sorted(
            pca_terms,
            key=lambda x: int(x.split("_")[1])
            if "_" in x and x.split("_")[1].isdigit()
            else x,
        )
        pca_terms = pca_terms[:4]

        # Build formula
        formula_covariates = [target_variable] + pca_terms
        formula_rhs = " + ".join(formula_covariates)
        formula = f"{weight_term} ~ {formula_rhs}"

        if debug:
            logger.debug("Formula for %s %s: %s", cell_type, modality, formula)

        try:
            model = smf.ols(formula=formula, data=covars_df)
            result = model.fit()

            # Create a dataframe with all coefficients for this cell type
            results_df = pd.DataFrame(
                {
                    "term": result.params.index,
                    "coefficient": result.params.values,
                    "stderr": result.bse.values,
                    "t-value": result.tvalues.values,
                    "p-value": result.pvalues.values,
                }
            )

            # Add metadata columns
            results_df.insert(0, "weight_term", weight_term)
            results_df.insert(0, "modality", modality)
            results_df.insert(0, "cell_type", cell_type)

            all_results.append(results_df)

        except Exception as e:
            logger.error(
                "Failed to run regression for %s %s: %s", cell_type, modality, str(e)
            )

    # Consolidate and save all results
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        out_file = (
            results_dir / f"{project}.all_{weight_term}_bias.{target_variable}.csv"
        )
        final_df.to_csv(out_file, index=False)
        logger.info("Successfully saved all regression results to %s", out_file)
    else:
        logger.warning("No results were generated.")


if __name__ == "__main__":
    main()
