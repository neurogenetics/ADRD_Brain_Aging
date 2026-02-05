from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet, get_dummies
import statsmodels.formula.api as smf
import concurrent.futures
import seaborn as sns
import matplotlib.pyplot as plt
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
    scale_dataframe,
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
        "--top-var-fraction",
        type=float,
        default=0.15,
        help="Fraction of top variable features to analyze (default: 0.15).",
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


def load_autosomal_features(features_file: Path, debug: bool = False) -> list[str]:
    logger.info(f"Loading features from {features_file}")
    features_df = read_csv(features_file)
    # Filter for autosomes (chr1-chr22)
    autosomes = [f"chr{i}" for i in range(1, 23)]
    autosomal_df = features_df[features_df["chr"].isin(autosomes)]

    if debug:
        logger.debug(f"Autosomal features shape: {autosomal_df.shape}")

    return autosomal_df["gene"].tolist()


def generate_latent_features(
    quants_df: DataFrame,
    covars_df: DataFrame,
    project: str,
    quants_dir: Path,
    out_figure_path: Path,
    title_suffix: str,
    top_var_fraction: float,
    imputer_type: str,
    debug: bool,
) -> DataFrame:
    logger.info("Begin modeling non-target variance in the data")

    # Load features and filter for autosomes
    features_file = quants_dir / f"{project}.features.csv"
    if features_file.exists():
        autosomal_genes = load_autosomal_features(features_file, debug)
        # Ensure we only use features present in the data
        # Note: quants_df features are in columns, samples in index
        candidate_features = quants_df.columns.intersection(autosomal_genes).tolist()
        logger.info(
            f"Restricted to {len(candidate_features)} autosomal features present in data"
        )
    else:
        logger.warning(
            f"Features file not found at {features_file}. Using all features."
        )
        candidate_features = quants_df.columns.tolist()

    variance_features = get_high_variance_features(
        quants_df[candidate_features], top_var_fraction
    )
    logger.info(f"Found {len(variance_features)} high variance features")
    max_count = int(
        min(
            quants_df[variance_features].shape[0], quants_df[variance_features].shape[1]
        )
        / 2
    )
    logger.info(f"Max components count: {max_count}")

    # deal with any missing values
    imputed_df = impute_missing_values(quants_df, variance_features, imputer_type)
    logger.info(f"Imputed DataFrame shape: {imputed_df.shape}")

    # Regress out age effects before PCA
    imputed_df = perform_regression_correction(imputed_df, covars_df, ["age"], debug)
    logger.info(
        f"Residuals DataFrame shape after regressing out age: {imputed_df.shape}"
    )

    pca_df = determine_pca_components(
        imputed_df, max_count, str(out_figure_path), debug, title_suffix
    )
    
    return pca_df


def fit_and_report_correlation(
    data_df: DataFrame,
    formula: str,
    description: str,
    return_tvalues: bool = False,
):
    """
    Helper to fit a GLM and report results.
    """
    logger.info(f"--- {description}: {formula} ---")
    try:
        model = smf.glm(formula=formula, data=data_df)
        result = model.fit()
        logger.info(result.summary())
        if return_tvalues:
            return result.tvalues
    except Exception as e:
        logger.warning(f"Failed to check correlations: {e}")
        return None


def check_correlations(
    data_df: DataFrame,
    target_var: str,
    covariates: list[str],
    description: str = "check correlations",
):
    """
    Checks if a target variable is correlated with a list of covariates.
    """
    covar_term_formula = " + ".join(covariates)
    this_formula = f"{target_var} ~ {covar_term_formula}"
    fit_and_report_correlation(data_df, this_formula, description)


def check_pcas_against_known_covariates(
    data_df: DataFrame,
    pca_cols: list[str],
    covariate_cols: list[str],
    out_prefix: Path = None,
    plot_title_suffix: str = "",
):
    covar_formula = " + ".join(covariate_cols)
    z_scores_list = []

    for pca in pca_cols:
        this_formula = f"{pca} ~ {covar_formula}"
        tvals = fit_and_report_correlation(
            data_df,
            this_formula,
            f"check if {pca} correlated with known covariates",
            return_tvalues=True,
        )

        if tvals is not None:
            tvals.name = pca
            z_scores_list.append(tvals)

    if out_prefix and z_scores_list:
        try:
            # Create DataFrame: rows = PCAs, cols = Covariates
            z_df = DataFrame(z_scores_list)

            # Drop Intercept if present
            if "Intercept" in z_df.columns:
                z_df = z_df.drop(columns=["Intercept"])

            # Create the heatmap
            plt.figure(figsize=(12, 10))
            # Use a diverging colormap, centering at 0
            sns.heatmap(z_df, cmap="RdBu_r", center=0, annot=True, fmt=".2f")
            plt.title(
                f"Z-statistics: PCA Components vs Known Covariates\n{plot_title_suffix}"
            )
            plt.xlabel("Known Covariates")
            plt.ylabel("PCA Components")
            plt.tight_layout()

            out_file = f"{out_prefix}_pca_covar_heatmap.png"
            plt.savefig(out_file)
            plt.close()
            logger.info(f"Saved PCA-Covariate heatmap to {out_file}")

        except Exception as e:
            logger.warning(f"Failed to generate heatmap: {e}")


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


def process_variance_results(
    results: dict,
    out_prefix: Path = None,
    suffix: str = "",
    debug: bool = False,
    plot_title_suffix: str = "",
):
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

        if out_prefix:
            try:
                # Prepare data for plotting
                plot_df = variance_fractions_df.melt(
                    var_name="Component", value_name="Variance Fraction"
                )

                plt.figure(figsize=(12, 6))
                sns.boxenplot(data=plot_df, x="Component", y="Variance Fraction")
                plt.xticks(rotation=45, ha="right")
                plt.title(
                    f"Variance Partitioning Results {suffix}\n{plot_title_suffix}"
                )
                plt.tight_layout()

                out_file = f"{out_prefix}_variance_boxen{suffix}.png"
                plt.savefig(out_file)
                plt.close()
                logger.info(f"Saved variance boxen plot to {out_file}")
            except Exception as e:
                logger.warning(f"Failed to create variance plot: {e}")
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


def perform_regression_correction(
    feature_df: DataFrame,
    covariate_df: DataFrame,
    covariate_cols: list[str],
    debug: bool = False,
) -> DataFrame:
    """
    Regresses out specified covariates from features.
    Handles categorical covariates via one-hot encoding.
    Imputes missing values in Y for fitting, but returns residuals based on original Y (preserving NaNs).
    """
    logger.info(f"Regressing out {covariate_cols} from features...")
    
    # Align indices
    common_idx = feature_df.index.intersection(covariate_df.index)
    if len(common_idx) < len(feature_df):
        logger.warning(
            f"Regression: Dropping {len(feature_df) - len(common_idx)} samples not in covariates."
        )

    Y_orig = feature_df.loc[common_idx]
    # Ensure covariates are a DataFrame
    X_source = covariate_df.loc[common_idx, covariate_cols]

    # One-hot encode covariates
    X = get_dummies(X_source, drop_first=True, dtype=float)
    if debug:
        peek_dataframe(X, "Encoded covariates matrix")

    # Handle missing in X (impute with mean)
    if X.isnull().values.any():
        logger.warning("Found missing values in covariates. Imputing with mean.")
        X = X.fillna(X.mean())

    # Handle missing in Y for fit (impute with mean)
    Y_fit = Y_orig.copy()
    if Y_fit.isnull().values.any():
        logger.warning("Found missing values in features. Imputing with mean for fit.")
        Y_fit = Y_fit.fillna(Y_fit.mean())

    # Fit
    reg = LinearRegression()
    reg.fit(X, Y_fit)

    # Residuals = Original - Predicted
    residuals = Y_orig - reg.predict(X)

    return residuals


def determine_pca_components(
    imputed_df: DataFrame,
    max_count: int,
    out_prefix: str = None,
    debug: bool = False,
    title_suffix: str = "",
) -> DataFrame:
    logger.info("Determine the number of PCA components to use")
    r2_values, rmse_values = iterate_model_component_counts(
        max_count, imputed_df, "PCA"
    )
    if debug:
        logger.debug(f"{r2_values=}")
        logger.debug(f"{rmse_values=}")

    # use max curvature of accuracy to select the number of components
    knee_rmse = component_from_max_curve(rmse_values, "RMSE", out_prefix, title_suffix)
    knee_r2 = component_from_max_curve(r2_values, "R2", out_prefix, title_suffix)
    num_comp = max(knee_rmse, knee_r2)
    logger.info(f"N = {num_comp} components will be used")

    pca_mdl, pca_df, _, _ = generate_selected_model(num_comp, imputed_df, "PCA")
    logger.info(f"PCA DataFrame shape: {pca_df.shape}")
    peek_dataframe(pca_df, "PCA variance components generated", debug)
    logger.info(f"Explained variance ratio: {pca_mdl.explained_variance_ratio_}")
    return pca_df


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    figs_dir = work_dir / "figures"
    
    # Ensure directories exist
    figs_dir.mkdir(parents=True, exist_ok=True)

    modality = args.modality
    cell_type = args.cell_type
    counts_term = f"{cell_type}_counts"
    probs_term = f"{cell_type}_probs" if modality == "atac" else None

    out_figure_path = figs_dir / f"{args.project}_{cell_type}_{modality}"
    title_suffix = f"{cell_type} ({modality.upper()})"

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

    # Define covariates for checking
    check_covariates = [
        "sex", "ancestry", "pmi", "ph", "smoker", "bmi", "pool", counts_term
    ]
    if modality == "atac":
        check_covariates.append(probs_term)

    check_correlations(covars_df, "age", check_covariates)

    # perform variance partition of known covariates
    fixed_effects = ["age", "bmi", "pmi", "ph", counts_term]
    if modality == "atac":
        fixed_effects.append(probs_term)

    random_effects = ["sex", "pool", "ancestry", "smoker"]

    results = run_variance_partition(
        data_df, quants_df.columns.tolist(), fixed_effects, random_effects, debug
    )

    process_variance_results(results, out_figure_path, "_known", debug, title_suffix)

    # Generate latent features representing non-target variance base on high variance features
    pca_df = generate_latent_features(
        quants_df,
        covars_df,
        args.project,
        quants_dir,
        out_figure_path,
        title_suffix,
        args.top_var_fraction,
        args.imputer_type,
        debug
    )

    # now redo perform variance partition of known covariates plus PCA components
    # extend the fixed effect terms to include the PCAs
    fixed_effects.extend(pca_df.columns.tolist())
    ext_data_df = data_df.merge(pca_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(ext_data_df, "Extended Data DataFrame", debug)

    check_correlations(
        ext_data_df, "age", pca_df.columns.tolist(), "check age vs PCA correlations"
    )

    known_covariates = [
        x for x in fixed_effects if x not in pca_df.columns
    ] + random_effects
    
    check_pcas_against_known_covariates(
        ext_data_df,
        pca_df.columns.tolist(),
        known_covariates,
        out_figure_path,
        title_suffix,
    )

    results = run_variance_partition(
        ext_data_df, quants_df.columns.tolist(), fixed_effects, random_effects, debug
    )

    process_variance_results(results, out_figure_path, "_final", debug, title_suffix)

    # Save final covariates terms to a file
    final_covariates = known_covariates + pca_df.columns.tolist()
    final_covariates_file = (
        info_dir / f"{args.project}.{cell_type}.{modality}.final_covariates.csv"
    )
    logger.info(f"Saving final covariates terms to {final_covariates_file}")
    ext_data_df[final_covariates].to_csv(final_covariates_file)

    # Regress out non-target covariates (everything except age)
    non_target_covariates = [x for x in final_covariates if x != "age"]
    
    # Features are the columns from the original quantified data
    feature_cols = quants_df.columns.tolist()

    residuals_df = perform_regression_correction(
        ext_data_df[feature_cols], ext_data_df, non_target_covariates, debug
    )

    # scale the dataframe features
    residuals_df = scale_dataframe(residuals_df)
    
    residuals_file = (
        quants_dir / f"{args.project}.{cell_type}.{modality}.residuals.parquet"
    )
    logger.info(f"Saving residuals to {residuals_file}")
    residuals_df.to_parquet(residuals_file)


if __name__ == "__main__":
    main()
