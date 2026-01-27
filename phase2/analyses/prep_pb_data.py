from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet, Series
import statsmodels.api as sm
import statsmodels.formula.api as smf
import concurrent.futures
import warnings

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


def perform_variance_partition(
    data_df: DataFrame,
    feature_name: str,
    fixed_effects: list[str],
    random_effects: list[str],
):
    """
    Mimics variancePartition by modeling continuous variables as fixed effects
    and categorical variables as random effects using statsmodels MixedLM.
    """
    try:
        # Create model columns list
        model_cols = fixed_effects + random_effects + [feature_name]

        # Create a temporary dataframe for modeling
        fit_df = data_df[model_cols].dropna()
        # Check if enough data remains
        if fit_df.empty:
            return None

        fit_df["dummy_group"] = 1

        # formula for fixed effects
        # Patsy handles basic formula parsing
        fixed_formula = f"Q('{feature_name}') ~ " + " + ".join(fixed_effects)

        # formula for random effects (categorical)
        vc_formula = {effect: f"0 + C({effect})" for effect in random_effects}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Fit the Linear Mixed Model
            model = smf.mixedlm(
                fixed_formula,
                data=fit_df,
                groups="dummy_group",
                re_formula="0",  # No random intercept for the dummy group itself
                vc_formula=vc_formula,
            )
            result = model.fit()
            # logger.info(result.summary()) # Silenced for batch processing

        # --- CALCULATE VARIANCE FRACTIONS ---

        # 1. Random Effects Variance (from vc_formula)
        try:
            vc_names = result.model.exog_vc.names
        except AttributeError:
            vc_names = [f"Var_Comp_{i}" for i in range(len(result.vcomp))]

        random_vars = dict(zip(vc_names, result.vcomp))

        # 2. Fixed Effects Variance
        # Var_explained = beta^2 * Var(X)
        fixed_vars = {}
        for term in fixed_effects:
            # Statsmodels params keys might differ slightly (e.g., Q('feature'))
            # But the predictors are in fixed_effects
            if term in result.params:
                beta = result.params[term]
                var_x = fit_df[term].var()
                fixed_vars[term] = (beta**2) * var_x

        # 3. Residual Variance
        residual_var = result.scale

        # 4. Combine and Normalize
        all_variances = {**random_vars, **fixed_vars, "Residual": residual_var}
        all_variances_series = Series(all_variances)

        total_var = all_variances_series.sum()
        fractions = all_variances_series / total_var

        return (feature_name, fractions)

    except Exception as e:
        logger.warning(f"Failed to process {feature_name}: {e}")
        return None


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
    fixed_effects = ["age", "bmi", "pmi", "ph", counts_term]
    if modality == "atac":
        fixed_effects.append(probs_term)

    random_effects = ["sex", "pool", "ancestry", "smoker"]

    # Identify all features (genes) from quants_df
    # Assuming columns in quants_df are the features
    features = quants_df.columns.tolist()
    logger.info(f"Starting variance partition for {len(features)} features...")

    results = {}

    # Use ProcessPoolExecutor for parallel processing
    # Using 'fork' context (default on Linux) makes passing data_df efficient (COW)
    with concurrent.futures.ProcessPoolExecutor() as executor:
        # Submit all tasks
        future_to_feature = {
            executor.submit(
                perform_variance_partition,
                data_df,
                feature,
                fixed_effects,
                random_effects,
            ): feature
            for feature in features
        }

        # Process results as they complete
        for i, future in enumerate(concurrent.futures.as_completed(future_to_feature)):
            feature = future_to_feature[future]
            try:
                res = future.result()
                if res is not None:
                    feat_name, fractions = res
                    results[feat_name] = fractions
            except Exception as exc:
                logger.debug(f"{feature} generated an exception: {exc}")

            if (i + 1) % 100 == 0:
                print(f"Processed {i + 1}/{len(features)}...", end="\r")

    print("\nProcessing complete.")

    if results:
        # Create DataFrame from results dictionary
        # Dictionary keys are index (genes), values are Series (columns)
        variance_fractions_df = DataFrame.from_dict(results, orient="index")

        print("\n--- Variance Fractions (Head) ---")
        peek_dataframe(
            variance_fractions_df, "computed variance partition fractions", debug
        )

        # Optionally describe the overall distribution
        print("\n--- Summary of Variance Explained ---")
        print(variance_fractions_df.describe())
    else:
        logger.warning("No results were generated.")


if __name__ == "__main__":
    main()
