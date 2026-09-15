import argparse
import logging
from pathlib import Path
import sys
import pandas as pd
import statsmodels.formula.api as smf

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Model the correlation between cell counts and target variable across all cell types to check for WLS weight bias."
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
        "--target-variable",
        type=str,
        default="age",
        help="The primary disease target variable column in covariates (default: 'age').",
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
    log_filename = f"{logs_dir}/{args.project}_batch_cell_counts_regression_{args.target_variable}.log"
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Logging configured. Writing to {log_filename}")

    results_dir.mkdir(parents=True, exist_ok=True)

    project = args.project
    weight_term = args.weight_term
    target_variable = args.target_variable

    # Find all final_covariates files for the project
    search_pattern = f"{project}.*.final_covariates.csv"
    covariates_files = list(info_dir.glob(search_pattern))
    
    if not covariates_files:
        logger.error(f"No final_covariates.csv files found in {info_dir} matching {search_pattern}")
        sys.exit(1)

    logger.info(f"Found {len(covariates_files)} covariate files to process.")

    # Pass 1: Compute total cell counts per sample per modality
    logger.info("Computing total cell counts per sample per modality...")
    modality_totals = {}
    for cov_file in covariates_files:
        parts = cov_file.name.split('.')
        if len(parts) < 5:
            continue
        modality = parts[-3]
        try:
            covars_df = pd.read_csv(cov_file, index_col=0)
            if weight_term in covars_df.columns:
                if modality not in modality_totals:
                    modality_totals[modality] = covars_df[weight_term].copy()
                else:
                    modality_totals[modality] = modality_totals[modality].add(covars_df[weight_term], fill_value=0)
        except Exception as e:
            logger.error(f"Failed to read {cov_file.name} for totals: {e}")

    all_results = []

    for cov_file in covariates_files:
        # Filename format: project.cell_type.modality.final_covariates.csv
        parts = cov_file.name.split('.')
        if len(parts) < 5:
            logger.warning(f"Skipping {cov_file.name}: does not match expected naming convention.")
            continue
        
        modality = parts[-3]
        cell_type = ".".join(parts[1:-3])

        logger.info(f"Processing cell_type: {cell_type}, modality: {modality}")

        covars_df = pd.read_csv(cov_file, index_col=0)
        
        if weight_term not in covars_df.columns:
            logger.warning(f"Weight term '{weight_term}' not found in {cov_file.name}. Skipping.")
            continue
            
        if target_variable not in covars_df.columns:
            logger.warning(f"Target variable '{target_variable}' not found in {cov_file.name}. Skipping.")
            continue

        # Compute percent cell counts
        percent_term = f"percent_{weight_term}"
        total_counts = modality_totals[modality].loc[covars_df.index]
        covars_df[percent_term] = (covars_df[weight_term] / total_counts) * 100

        # Identify PCA terms, limiting to the first 4 to match target regression modeling
        pca_terms = [col for col in covars_df.columns if col.startswith("PCA_")]
        
        # Sort to ensure stable order (e.g., PCA_0, PCA_1, PCA_2, PCA_3)
        pca_terms = sorted(pca_terms, key=lambda x: int(x.split('_')[1]) if '_' in x and x.split('_')[1].isdigit() else x)
        pca_terms = pca_terms[:4]
        
        formula_covariates = [target_variable] + pca_terms
        formula_rhs = " + ".join(formula_covariates)
        
        # Run regression for both raw counts and percent counts
        for term_to_test in [weight_term, percent_term]:
            formula = f"{term_to_test} ~ {formula_rhs}"
            
            if debug:
                logger.debug(f"Formula for {cell_type} {modality}: {formula}")
            
            try:
                model = smf.ols(formula=formula, data=covars_df)
                result = model.fit()
                
                # Create a dataframe with all coefficients for this cell type
                results_df = pd.DataFrame({
                    "term": result.params.index,
                    "coefficient": result.params.values,
                    "stderr": result.bse.values,
                    "t-value": result.tvalues.values,
                    "p-value": result.pvalues.values
                })
                
                # Add metadata columns
                results_df.insert(0, "weight_term", term_to_test)
                results_df.insert(0, "modality", modality)
                results_df.insert(0, "cell_type", cell_type)
                
                all_results.append(results_df)
                
            except Exception as e:
                logger.error(f"Failed to run regression for {cell_type} {modality} ({term_to_test}): {e}")

    # Consolidate and save all results
    if all_results:
        final_df = pd.concat(all_results, ignore_index=True)
        out_file = results_dir / f"{project}.all_{weight_term}_bias.{target_variable}.csv"
        final_df.to_csv(out_file, index=False)
        logger.info(f"Successfully saved all regression results to {out_file}")
    else:
        logger.warning("No results were generated.")

if __name__ == "__main__":
    main()
