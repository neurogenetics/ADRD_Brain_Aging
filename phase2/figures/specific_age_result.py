import argparse
import logging
from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scanpy import read_h5ad

# Suppress warnings
warnings.simplefilter('ignore')

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run specific age result analysis and generate figures."
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
        default="rna",
        choices=["rna", "atac"],
        help="Data modality (rna or atac).",
    )
    parser.add_argument(
        "--feature",
        type=str,
        required=True,
        help="Feature name (gene or peak) to analyze.",
    )
    parser.add_argument(
        "--cell-type",
        type=str,
        required=True,
        help="Cell type to analyze.",
    )
    parser.add_argument(
        "--regression-type",
        type=str,
        default="glm_tweedie",
        choices=["glm", "glm_tweedie", "rlm"],
        help="Regression method to use.",
    )
    parser.add_argument(
        "--debug", 
        action="store_true", 
        help="Enable debug output."
    )
    return parser.parse_args()


def main():
    args = parse_args()
    debug = args.debug

    # Configure logging
    log_filename = f"{args.modality}_{args.cell_type}_{args.feature}_specific_analysis.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ],
        force=True
    )
    logger.info(f"Logging configured. Writing to {log_filename}")

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    results_dir = work_dir / "results"
    figures_dir = work_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    project = args.project
    modality = args.modality.lower() # ensure lowercase to match file patterns in scripts
    feature = args.feature
    cell_type = args.cell_type
    regression_type = args.regression_type

    # Note: Notebook logic for 'prefix_type' (broad vs specific) was based on 'category'
    # 'category' was 'curated_type' (broad) or 'cluster_name' (specific).
    # Since we are automating, we might need to know which one it is.
    # However, existing scripts (like post_pseudobulk) assumed everything is consistent.
    # We'll try to guess or use a default. The notebook default was 'curated_type' -> 'broad'.
    # If the user provides a cell_type that is specific, this might break if we hardcode 'broad'.
    # Let's check for the file existence to decide if possible, or just default to what notebook had.
    # The notebook had: category = 'curated_type' -> prefix_type = 'broad'.
    # Let's try 'broad' first, if file not found, maybe try 'specific'? 
    # Or better, just rely on the file path constructed.
    # We will assume 'broad' (curated_type) as per notebook default.
    prefix_type = 'broad' 
    
    # 1. Read Regression Results
    # File format: {project}.{modality}.{prefix_type}.{cell_type}.{regression_type}.age.csv
    # NOTE: Previous script post_pseudobulk used: {project}.{modality}.all_celltypes.{regression_type}.age.csv
    # If this script is meant to look at the aggregated results, we should use that.
    # But the notebook looked at individual files: 
    # in_file = f'{results_dir}/{project}.{modality}.{prefix_type}.{cell_type}.{REGRESSION_TYPE}.age.csv'
    # We will try to read from the aggregated 'all_celltypes' file first as it's more robust now.
    
    all_results_file = results_dir / f'{project}.{modality}.all_celltypes.{regression_type}.age.csv'
    if all_results_file.exists():
        logger.info(f"Reading aggregated results from {all_results_file}")
        glm_results = pd.read_csv(all_results_file)
        # Filter for specific feature and cell type
        specific_result = glm_results.loc[
            (glm_results.feature == feature) & 
            (glm_results.tissue == cell_type)
        ]
        if not specific_result.empty:
             logger.info(f"Found result for {feature} in {cell_type}:")
             print(specific_result)
        else:
            logger.warning(f"No result found for {feature} in {cell_type} in aggregated file.")
    else:
        logger.warning(f"Aggregated results file not found: {all_results_file}")

    # 2. Load Quantified Data (Pseudobulk)
    # File: {quants_dir}/{project}.{cell_type}.{modality}.parquet
    pb_file = quants_dir / f'{project}.{cell_type}.{modality}.parquet'
    logger.info(f"Reading pseudobulk data from {pb_file}")
    
    if not pb_file.exists():
        logger.error(f"Pseudobulk file not found: {pb_file}")
        return

    quants_df = pd.read_parquet(pb_file)
    logger.info(f"Quantifications shape: {quants_df.shape}")
    
    if feature not in quants_df.columns:
        logger.error(f"Feature {feature} not found in quantifications.")
        return

    # 3. Load Covariates (from CSV)
    # File: sample_info/{project}.{cell_type}.{modality}.final_covariates.csv
    covars_file = work_dir / 'sample_info' / f'{project}.{cell_type}.{modality}.final_covariates.csv'
    logger.info(f"Reading covariates from {covars_file}")
    
    if not covars_file.exists():
        logger.error(f"Covariates file not found: {covars_file}")
        return
        
    covars_df = pd.read_csv(covars_file, index_col=0)
    logger.info(f"Covariates shape: {covars_df.shape}")

    # 4. Prepare Dataframe
    # Check for 'cell_counts' and rename to 'cell_count'
    if 'cell_counts' in covars_df.columns:
        covars_df = covars_df.rename(columns={'cell_counts': 'cell_count'})
    
    # Merge
    # quants_df should be just features
    # covars_df has all covariates including cell_count
    
    data_df = quants_df[[feature]].merge(covars_df, how='inner', left_index=True, right_index=True)
    logger.info(f"Merged data shape: {data_df.shape}")
    
    # 5. Run Regression (Re-verification)
    covariate_terms = ['sex', 'ancestry', 'pmi', 'ph', 'smoker', 'bmi', 'pool']
    # Filter terms present in data_df
    covariate_terms = [t for t in covariate_terms if t in data_df.columns]
    covar_term_formula = ' + '.join(covariate_terms)
    
    endo_term = feature
    exog_term = 'age'
    
    # Check for Tweedie vs GLM vs RLM
    # Formula: feature ~ age + covariates + cell_count
    # Note: GLM Tweedie in notebook used Q("feature") ~ ...
    
    formula = f'Q("{endo_term}") ~ {exog_term} + {covar_term_formula} + cell_count'
    
    logger.info(f"Running {regression_type} for {feature}...")
    
    try:
        if regression_type == 'glm_tweedie':
            model = smf.glm(
                formula=formula, 
                data=data_df, 
                family=sm.families.Tweedie(link=sm.families.links.log(), var_power=1.6, eql=True)
            )
        elif regression_type == 'glm':
            model = smf.glm(formula=formula, data=data_df)
        elif regression_type == 'rlm':
            model = smf.rlm(formula=formula, data=data_df)
        else:
            logger.error(f"Unknown regression type: {regression_type}")
            return

        result = model.fit()
        print(result.summary())
    except Exception as e:
        logger.error(f"Regression failed: {e}")

    # 6. Generate Figures
    logger.info("Generating figures...")
    
    # Distribution Plot
    plt.figure(figsize=(8, 6))
    sns.histplot(data_df[feature], kde=True)
    plt.title(f"Distribution of {feature} in {cell_type}")
    plt.savefig(figures_dir / f"{project}.{modality}.{cell_type}.{feature}.dist.png")
    plt.close()
    
    # Scatter Plot: Age vs Feature (Size by Cell Count)
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='age', y=feature, size='cell_count', data=data_df)
    plt.title(f"{feature} vs Age in {cell_type}")
    plt.savefig(figures_dir / f"{project}.{modality}.{cell_type}.{feature}.scatter.png")
    plt.close()

    # Regression Plot (Robust)
    plt.figure(figsize=(8, 6))
    sns.lmplot(x='age', y=feature, data=data_df, robust=True, height=6, aspect=1.2)
    plt.title(f"{feature} vs Age (Robust Fit)")
    plt.savefig(figures_dir / f"{project}.{modality}.{cell_type}.{feature}.lm_robust.png")
    plt.close()
    
    # Regression Plot (Standard)
    plt.figure(figsize=(8, 6))
    sns.lmplot(x='age', y=feature, data=data_df, height=6, aspect=1.2)
    plt.title(f"{feature} vs Age (Standard Fit)")
    plt.savefig(figures_dir / f"{project}.{modality}.{cell_type}.{feature}.lm.png")
    plt.close()
    
    logger.info("Done.")


if __name__ == "__main__":
    main()
