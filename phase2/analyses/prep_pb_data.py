from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet, Series
import statsmodels.api as sm
import statsmodels.formula.api as smf

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert single-cell AnnData to pseudobulk profiles."
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
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def peek_dataframe(df: DataFrame, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    print(f"{df.shape=}")
    if verbose:
        if len(df.columns) < 25:
            print(tabulate(df.head(), headers="keys", tablefmt="psql"))
        else:
            print(f"{df.index.values[0:15]=}")
            print(f"{df.columns.values[0:15]=}")


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"

    modality = "rna"
    cell_type = "Microglia"
    counts_term = f"{cell_type}_counts"
    if modality == "atac":
        probs_term = f"{cell_type}_probs"

    # read the covariate file
    covars_file = info_dir / f"{args.project}.covariates.{modality}.csv"
    covars_df = read_csv(covars_file, index_col=0)
    peek_dataframe(covars_df, f"loaded the covariates file: {covars_file}", debug)
    if debug:
        print(covars_df.sex.value_counts())
        print(covars_df.ancestry.value_counts())
        print(covars_df.pool.value_counts())
        print(covars_df.pmi.describe())
        print(covars_df.ph.describe())
        print(covars_df.smoker.describe())
        print(covars_df.bmi.describe())
        print(covars_df[counts_term].describe())
        if modality == "atac":
            print(covars_df[probs_term].describe())
        print(covars_df.info())

    # load the quantified data
    data_file = quants_dir / f"{args.project}.{cell_type}.{modality}.parquet"
    quants_df = read_parquet(data_file)
    peek_dataframe(quants_df, f"loaded the quantified data file: {data_file}", debug)

    # check if age, exogenous variable, is correlated with any ouf the covariate terms
    covariate_terms = [
        "sex",
        "ancestry",
        "pmi",
        "ph",
        "smoker",
        "bmi",
        "pool",
        counts_term,
    ]
    if modality == "atac":
        covariate_terms.append(probs_term)
    covar_term_formula = " + ".join(covariate_terms)
    this_formula = f"age ~ {covar_term_formula}"
    logger.info(
        f"--- check if age correlated with other covariates: {covar_term_formula} ---"
    )
    model = smf.glm(formula=this_formula, data=covars_df)
    result = model.fit()
    logger.info(result.summary())

    # combine the quantifications and covariates for modeling and cleaning of non-target variance
    data_df = covars_df.merge(quants_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(data_df, "merged the covariates and quantifications", debug)

    # perform variance partition of known covariates
    # Create a temporary dataframe for modeling to avoid modifying the main data_df
    # and to ensure statsmodels doesn't fail due to NaN alignment issues.
    # model_cols = ["CHURC1", "sex", "pool", "pmi", "ph"]
    model_cols = covariate_terms + ["CHURC1"]
    fit_df = data_df[model_cols].dropna()
    fit_df["dummy_group"] = 1

    # Fit the Linear Mixed Model
    # Continuous variables are modeled as Fixed Effects (main formula).
    # Categorical variables are modeled as Random Effects (vc_formula).
    # This mimics variancePartition behavior more accurately.
    
    # Continuous variables for fixed effects
    fixed_effects = ["pmi", "ph", counts_term]
    # "ancestry" and "smoker" are categorical but sometimes treated as fixed if levels are few; 
    # however, variancePartition often treats batches/individuals as random.
    # Based on your previous code, ancestry/smoker were in vc_formula (random).
    # If they are categorical, keep them in vc_formula. If they are numeric/continuous, move here.
    # Assuming standard covariates:
    # - pmi, ph, counts: Continuous -> Fixed
    # - sex, pool, ancestry, smoker: Categorical -> Random
    
    fixed_formula = "CHURC1 ~ " + " + ".join(fixed_effects)
    
    # Categorical variables for random effects
    vc_formula = {
        "sex": "0 + C(sex)",
        "pool": "0 + C(pool)",
        "ancestry": "0 + C(ancestry)",
        "smoker": "0 + C(smoker)",
    }

    model = smf.mixedlm(
        fixed_formula,
        data=fit_df,
        groups="dummy_group",
        re_formula="0", # No random intercept for the dummy group itself
        vc_formula=vc_formula,
    )
    result = model.fit()
    logger.info(result.summary())

    # --- CALCULATE VARIANCE FRACTIONS ---
    
    # 1. Random Effects Variance (from vc_formula)
    # result.vcomp contains the variance estimates for the random terms as an array
    # We need to map these to their names
    try:
        vc_names = result.model.exog_vc.names
    except AttributeError:
        vc_names = [f"Var_Comp_{i}" for i in range(len(result.vcomp))]
        
    random_vars = dict(zip(vc_names, result.vcomp))
    
    # 2. Fixed Effects Variance
    # For each fixed effect X with coefficient beta, Var_explained = Var(beta * X)
    # Since beta is a constant, Var(beta * X) = beta^2 * Var(X)
    fixed_vars = {}
    for term in fixed_effects:
        # Get coefficient (params) for the term
        if term in result.params:
            beta = result.params[term]
            # Calculate variance of the predictor in the data
            var_x = fit_df[term].var()
            fixed_vars[term] = (beta ** 2) * var_x
            
    # 3. Residual Variance
    residual_var = result.scale
    
    # 4. Combine and Normalize
    all_variances = {**random_vars, **fixed_vars, "Residual": residual_var}
    all_variances_series = Series(all_variances)
    
    total_var = all_variances_series.sum()
    fractions = all_variances_series / total_var
    
    print("\n--- Variance Fractions ---")
    print(fractions.sort_values(ascending=False))


if __name__ == "__main__":
    main()
