from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet, Series
import statsmodels.api as sm
import statsmodels.formula.api as smf
import concurrent.futures
import warnings
from variance_utils import (
    get_high_variance_features,
    perform_variance_partition,
    iterate_model_component_counts,
    component_from_max_curve,
    generate_selected_model,
)

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
    parser.add_argument(
        "--top-var-percent",
        type=float,
        default=0.15,
        help="Percentage of top variable features to analyze (default: 0.15).",
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
        variance_fractions_df = DataFrame.from_dict(results, orient="index").round(4)

        print("\n--- Variance Fractions (Head) ---")
        peek_dataframe(
            variance_fractions_df, "computed variance partition fractions", debug
        )

        # Optionally describe the overall distribution
        print("\n--- Summary of Variance Explained ---")
        print(variance_fractions_df.describe())
    else:
        logger.warning("No results were generated.")

    # Generate latent features representing non-target variance base on high variance features
    logger.info("Begin modeling non-target variance in the data")
    variance_features = get_high_variance_features(quants_df)
    logger.info(f"found {len(variance_features)} high variance features")
    max_count = int(
        min(
            quants_df[variance_features].shape[0], quants_df[variance_features].shape[1]
        )
        / 2
    )
    print(f"max count is {max_count}")

    # deal with any missing values
    logger.info("impute missing values")
    from sklearn.experimental import enable_iterative_imputer
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer

    imputer_type = "knn"
    if imputer_type == "iterative":
        logger.info("imputing with IterativeImputer")
        imputer = IterativeImputer(
            estimator=ExtraTreesRegressor(n_estimators=10, random_state=42, n_jobs=-1),
            max_iter=20,
            random_state=42,
            n_nearest_features=20,
        )
    elif imputer_type == "knn":
        neighbor_cnt = int(len(quants_df) ** 0.5)
        logger.info(f"imputing with KNNImputer with n = {neighbor_cnt}")
        imputer = KNNImputer(n_neighbors=neighbor_cnt)
    else:  # simple
        logger.info("imputing with SimpleImputer")
        imputer = SimpleImputer(strategy="mean")
    imputer.set_output(transform="pandas")
    imputed_df = imputer.fit_transform(quants_df[variance_features])
    print(imputed_df)

    logger.info("determine the number of PCA components to use")
    r2_values, rmse_values = iterate_model_component_counts(
        max_count, imputed_df, "PCA"
    )
    if debug:
        print(f"{r2_values=}")
        print(f"{rmse_values=}")

    # use max curvature of accuracy to select the number of components
    knee_rmse = component_from_max_curve(rmse_values, "RMSE")
    knee_r2 = component_from_max_curve(r2_values, "R2")
    num_comp = max(knee_rmse, knee_r2)
    logger.info(f"N = {num_comp} components will be used")

    pca_mdl, pca_df, _, _ = generate_selected_model(num_comp, imputed_df, "PCA")
    print(f"shape of pca_df is {pca_df.shape}")
    peek_dataframe(pca_df, "PCA variance components generated", debug)
    # plot_pca_pair(pca_df.merge(covs_df, how='left',
    #                         left_index=True, right_index=True),
    #             'PCA_0', 'PCA_1', hue_cov='pool_num')
    print(pca_mdl.explained_variance_ratio_)


if __name__ == "__main__":
    main()
