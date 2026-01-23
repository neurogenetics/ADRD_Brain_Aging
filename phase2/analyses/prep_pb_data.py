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
        print(covars_df.info())

    # load the quantified data
    data_file = quants_dir / f"{args.project}.{cell_type}.{modality}.parquet"
    quants_df = read_parquet(data_file)
    peek_dataframe(quants_df, f"loaded the quantified data file: {data_file}", debug)

    # check if age, exogenous variable, is correlated with any ouf the covariate terms
    covariate_terms = ["sex", "ancestry", "pmi", "ph", "smoker", "bmi", "pool"]
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
    # Fit the Linear Mixed Model
    # "expression ~ 1" means we are modeling the intercept (mean)
    # "groups" defines the random effect (e.g., individual)
    # variancePartition fits multiple random effects; statsmodels handles this via 'vc_formula'
    model = smf.mixedlm(
        "CHURC1 ~ 1",
        data=data_df,
        groups=data_df.index,
        vc_formula={
            "sex": "0 + C(sex)",
            "pool": "0 + C(pool)",
            "pmi": "0 + pmi",
            "ph": "0 + ph",
            # "ancestry": "0 + C(ancestry)",
            # "smoker": "0 + smoker",
        },
    )
    result = model.fit()
    logger.info(result.summary())
    # 2. EXTRACT VARIANCE COMPONENTS (Robust method)
    # Get the raw variance components array
    var_components = result.vcomp
    # Try to get the names of these components from the model
    # (These usually correspond to the order in the array: Group var, then vc_formula vars)
    try:
        # This attribute exists in most statsmodels versions to label the VC params
        vc_names = result.model.exog_vc.names
    except AttributeError:
        # Fallback if names can't be found automatically
        vc_names = [f"Var_Comp_{i}" for i in range(len(var_components))]
    # 3. Create a labeled Pandas Series
    # This ensures 'fractions' is a Series, not a raw numpy array
    fractions = Series(var_components, index=vc_names)
    # 4. Calculate Fractions
    residual_var = result.scale
    total_var = fractions.sum() + residual_var
    fractions = fractions / total_var
    # Now this will work because 'fractions' is a pandas Series, not a numpy array
    fractions["Residual"] = residual_var / total_var
    print(fractions)


if __name__ == "__main__":
    main()
