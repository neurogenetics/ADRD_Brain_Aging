from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet
import statsmodels.formula.api as smf
import concurrent.futures
from sklearn.experimental import enable_iterative_imputer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.linear_model import LinearRegression
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
        "--modality",
        type=str,
        default="rna",
        help="Data modality (e.g., rna, atac).",
    )
    parser.add_argument(
        "--cell-type",
        type=str,
        default="Microglia",
        help="Cell type to analyze.",
    )
    parser.add_argument(
        "--imputer-type",
        type=str,
        default="knn",
        choices=["iterative", "knn", "simple"],
        help="Type of imputation to use.",
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
    logger.info(f"DataFrame shape: {df.shape}")
    if verbose:
        if len(df.columns) < 25:
            print(tabulate(df.head(), headers="keys", tablefmt="psql"))
        else:
            print(f"{df.index.values[0:15]=}")
            print(f"{df.columns.values[0:15]=}")


def load_covariates(
    info_dir: Path, project: str, modality: str, debug: bool = False
) -> DataFrame:
    covars_file = info_dir / f"{project}.covariates.{modality}.csv"
    covars_df = read_csv(covars_file, index_col=0)
    peek_dataframe(covars_df, f"Loaded covariates file: {covars_file}", debug)
    return covars_df


def load_quantified_data(
    quants_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> DataFrame:
    data_file = quants_dir / f"{project}.{cell_type}.{modality}.parquet"
    quants_df = read_parquet(data_file)
    peek_dataframe(quants_df, f"Loaded quantified data file: {data_file}", debug)
    return quants_df


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    figs_dir = work_dir / "figures"

    modality = args.modality
    cell_type = args.cell_type
    counts_term = f"{cell_type}_counts"
    probs_term = f"{cell_type}_probs" if modality == "atac" else None

    # Load data
    covars_df = load_covariates(info_dir, args.project, modality, debug)
    if debug:
        logger.debug(covars_df.describe())
        logger.debug(covars_df.info())

    quants_df = load_quantified_data(
        quants_dir, args.project, cell_type, modality, debug
    )

    # combine the quantifications and covariates for modeling and cleaning of non-target variance
    data_df = covars_df.merge(quants_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(data_df, "merged the covariates and quantifications", debug)

    check_covariate_correlations(covars_df, counts_term, probs_term, modality, "age")

    # perform variance partition of known covariates
    fixed_effects = ["age", "bmi", "pmi", "ph", counts_term]
    if modality == "atac":
        fixed_effects.append(probs_term)

    random_effects = ["sex", "pool", "ancestry", "smoker"]

    results = run_variance_partition(
        data_df, quants_df.columns.tolist(), fixed_effects, random_effects, debug
    )

    process_variance_results(results, debug)

    # Generate latent features representing non-target variance base on high variance features
    logger.info("Begin modeling non-target variance in the data")
    variance_features = get_high_variance_features(quants_df)
    logger.info(f"Found {len(variance_features)} high variance features")
    max_count = int(
        min(
            quants_df[variance_features].shape[0], quants_df[variance_features].shape[1]
        )
        / 2
    )
    logger.info(f"Max components count: {max_count}")

    # deal with any missing values
    imputed_df = impute_missing_values(quants_df, variance_features, args.imputer_type)
    logger.info(f"Imputed DataFrame shape: {imputed_df.shape}")

    # Regress out age effects before PCA
    imputed_df = regress_out_covariate(imputed_df, covars_df, "age")
    logger.info(
        f"Residuals DataFrame shape after regressing out age: {imputed_df.shape}"
    )

    out_figure_path = figs_dir / f"{args.project}_{cell_type}_{modality}"
    pca_df = determine_pca_components(imputed_df, max_count, out_figure_path, debug)

    # now redo perform variance partition of known covariates plus PCA components
    # extend the fixed effect terms to include the PCAs
    fixed_effects.extend(pca_df.columns.tolist())
    ext_data_df = data_df.merge(pca_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(ext_data_df, "Extended Data DataFrame", debug)

    results = run_variance_partition(
        ext_data_df, quants_df.columns.tolist(), fixed_effects, random_effects, debug
    )

    process_variance_results(results, debug)


def check_covariate_correlations(
    covars_df: DataFrame,
    counts_term: str,
    probs_term: str,
    modality: str,
    target_var: str = "age",
):
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
    this_formula = f"{target_var} ~ {covar_term_formula}"
    logger.info(
        f"--- check if {target_var} correlated with other covariates: {covar_term_formula} ---"
    )
    try:
        model = smf.glm(formula=this_formula, data=covars_df)
        result = model.fit()
        logger.info(result.summary())
    except Exception as e:
        logger.warning(f"Failed to check correlations: {e}")


def run_variance_partition(
    data_df: DataFrame,
    features: list,
    fixed_effects: list,
    random_effects: list,
    debug: bool = False,
) -> dict:
    logger.info(f"Starting variance partition for {len(features)} features...")
    results = {}

    # Use ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor() as executor:
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
    return results


def process_variance_results(results: dict, debug: bool = False):
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


def impute_missing_values(
    quants_df: DataFrame, variance_features: list, imputer_type: str
) -> DataFrame:
    logger.info("impute missing values")

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
    return imputed_df


def regress_out_covariate(
    data_df: DataFrame, covars_df: DataFrame, covariate: str
) -> DataFrame:
    logger.info(f"Regressing out {covariate} from features...")
    # intersect indices
    common_idx = data_df.index.intersection(covars_df.index)
    if len(common_idx) < len(data_df):
        logger.warning(
            f"Dropping {len(data_df) - len(common_idx)} samples missing {covariate} info."
        )

    y = data_df.loc[common_idx]
    x = covars_df.loc[common_idx, [covariate]]

    reg = LinearRegression()
    reg.fit(x, y)
    residuals = y - reg.predict(x)

    return residuals


def determine_pca_components(
    imputed_df: DataFrame, max_count: int, out_prefix: str = None, debug: bool = False
) -> DataFrame:
    logger.info("Determine the number of PCA components to use")
    r2_values, rmse_values = iterate_model_component_counts(
        max_count, imputed_df, "PCA"
    )
    if debug:
        logger.debug(f"{r2_values=}")
        logger.debug(f"{rmse_values=}")

    # use max curvature of accuracy to select the number of components
    knee_rmse = component_from_max_curve(rmse_values, "RMSE", out_prefix)
    knee_r2 = component_from_max_curve(r2_values, "R2", out_prefix)
    num_comp = max(knee_rmse, knee_r2)
    logger.info(f"N = {num_comp} components will be used")

    pca_mdl, pca_df, _, _ = generate_selected_model(num_comp, imputed_df, "PCA")
    logger.info(f"PCA DataFrame shape: {pca_df.shape}")
    peek_dataframe(pca_df, "PCA variance components generated", debug)
    logger.info(f"Explained variance ratio: {pca_mdl.explained_variance_ratio_}")
    return pca_df


if __name__ == "__main__":
    main()
