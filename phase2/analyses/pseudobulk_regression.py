import argparse
import logging
import random
import warnings
from multiprocessing import Process
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from pandas import DataFrame as PandasDF
from pandas import read_csv, read_parquet

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    force=True
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run pseudobulk regression analysis."
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


def load_final_covariates(
    info_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> PandasDF:
    final_covars_file = info_dir / f"{project}.{cell_type}.{modality}.final_covariates.csv"
    if not final_covars_file.exists():
        logger.error(f"Final covariates file not found: {final_covars_file}")
        raise FileNotFoundError(f"Final covariates file not found: {final_covars_file}")

    final_covars_df = read_csv(final_covars_file, index_col=0)
    if debug:
        logger.info(f"Loaded final covariates file: {final_covars_file}")
        print(final_covars_df.head())
    return final_covars_df


def load_residuals(
    quants_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> PandasDF:
    data_file = quants_dir / f"{project}.{cell_type}.{modality}.residuals.parquet"
    if not data_file.exists():
        logger.error(f"Residuals data file not found: {data_file}")
        raise FileNotFoundError(f"Residuals data file not found: {data_file}")

    quants_df = read_parquet(data_file)
    if debug:
        logger.info(f"Loaded residuals data file: {data_file}")
        print(quants_df.head())
    return quants_df


def glm_model(formula: str, df: PandasDF, model_type: str = 'rlm'):
    if model_type == 'glm_tweedie':
        model = smf.glm(formula=formula, data=df, 
                        family=sm.families.Tweedie(link=sm.families.links.Log(), 
                                                   var_power=1.6, 
                                                   eql=True))
    elif model_type == 'rlm':
        model = smf.rlm(formula=formula, data=df)        
    elif model_type == 'glm':
        model = smf.glm(formula=formula, data=df)        
    result = model.fit()
    return result


def glm_age(df: PandasDF, feature: str, regression_type: str, verbose: bool = False) -> list:
    endo_term = feature
    exog_term = 'age'
    # Note: Using only 'age' as the model term since other covariates were regressed out in prep_pb_data.py
    # However, we need to make sure 'age' is in the dataframe.
    # The residuals dataframe should have features as columns and samples as rows.
    # The 'age' variable comes from the merged covariates dataframe.
    
    this_formula = f'Q("{endo_term}") ~ {exog_term}'
    
    try:
        # run GLM via statsmodel
        result = glm_model(this_formula, df[[endo_term, exog_term]], model_type=regression_type)
        ret_list = [endo_term, result.params['Intercept'], 
                    result.params[exog_term], result.bse[exog_term], 
                    result.tvalues[exog_term], result.pvalues[exog_term]]
        if verbose:
            print(f'df shape {df.shape}')
            print(result.summary())
            print(['feature', 'intercept', 'coef', 'stderr', 'z', 'p-value'])
            print(ret_list)
    except Exception as e:
        if verbose:
            print(f'Caught Error for {endo_term}: {e}')
        ret_list = [endo_term] + [np.nan] * 5
  
    return ret_list


def regress_age(quants: PandasDF, covars: PandasDF, cell_name: str, 
                modality: str, regression_type: str, results_dir: Path, project: str) -> None:
    
    logger.info(f"Starting regression for {cell_name} using {regression_type}")
    
    # Merge residuals with covariates (we need 'age')
    # Note: covars contains 'age' and other technical factors.
    # Even though we regressed out technical factors, we still need 'age' for the final association test.
    data_df = quants.merge(covars[['age']], how='inner', left_index=True, right_index=True)
    
    features = quants.columns.tolist()
    logger.info(f"Processing {len(features)} features...")
    
    type_results = []
    for i, feature in enumerate(features):
        type_results.append(glm_age(data_df, feature, regression_type))
        if (i + 1) % 1000 == 0:
            print(f"Processed {i + 1}/{len(features)}...", end="\r")
            
    print("\nProcessing complete.")

    results_df = PandasDF(data=type_results,
                          columns=['feature', 'intercept', 'coef', 'stderr', 
                                    'z', 'p-value'])
    results_df['tissue'] = cell_name
    
    # Save results
    out_file = results_dir / f'{project}.{modality}.{cell_name}.{regression_type}.age.csv'
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved results to {out_file}")


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    results_dir = work_dir / "results"
    
    # Ensure results directory exists
    results_dir.mkdir(parents=True, exist_ok=True)

    modality = args.modality
    cell_type = args.cell_type
    project = args.project
    
    # Load final covariates (contains 'age' and we might need it, 
    # but strictly speaking we just need 'age' which is in the original covariates too.
    # Using final_covariates is safer to ensure sample alignment if any filtering happened).
    covars_df = load_final_covariates(info_dir, project, cell_type, modality, debug)
    
    # Load residuals (clean data)
    residuals_df = load_residuals(quants_dir, project, cell_type, modality, debug)

    # Run regression
    regress_age(
        residuals_df, 
        covars_df, 
        cell_type, 
        modality, 
        args.regression_type, 
        results_dir, 
        project
    )


if __name__ == "__main__":
    main()
