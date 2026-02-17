from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet, get_dummies
from sklearn.experimental import enable_iterative_imputer
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.linear_model import LinearRegression
import numpy as np
from prince import FAMD
from variance_utils import (
    get_high_variance_features,
    iterate_model_component_counts,
    component_from_max_curve,
    generate_selected_model,
    scale_dataframe,
    check_correlations,
)

# Configure logging
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
        default=0.25,
        help="Fraction of top variable features to analyze (default: 0.25).",
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
    covariates_df: DataFrame,
    covariate_cols: list[str],
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

    # Regress out known covariates effects before PCA
    imputed_df = perform_regression_correction(
        imputed_df, covariates_df, covariate_cols, debug
    )
    logger.info(
        f"Residuals DataFrame shape after regressing out known covariates: {imputed_df.shape}"
    )

    pca_df = determine_pca_components(
        imputed_df, max_count, str(out_figure_path), debug, title_suffix
    )

    return pca_df


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


def generate_famd_features(
    covariates_df: DataFrame,
    covariate_cols: list[str],
    out_prefix: str,
    title_suffix: str,
    n_components: int = 5,
    debug: bool = False,
) -> DataFrame:
    """
    Generate FAMD components from mixed continuous/categorical covariates.
    Uses a fixed number of components (default 5) for covariate compression.
    """
    logger.info(f"Generating {n_components} FAMD features from final covariates...")

    # Filter for specified covariates
    famd_input = covariates_df[covariate_cols].copy()

    if debug:
        peek_dataframe(famd_input, "FAMD Input Data")

    # Ensure n_components is valid given data dimensions
    max_possible = min(famd_input.shape[0], famd_input.shape[1])
    if n_components > max_possible:
        logger.warning(
            f"Requested {n_components} FAMD components, but data dimensions limit to {max_possible}. "
            f"Using {max_possible}."
        )
        n_components = max_possible

    try:
        # Fit final model
        logger.info(f"Running Prince FAMD with n_components={n_components}...")
        final_famd = FAMD(
            n_components=n_components,
            n_iter=3,  # Increased iterations for stability
            random_state=42,
            engine="sklearn",
            handle_unknown="error",
        )
        # FAMD return type can vary (DataFrame or array-like)
        famd_coords = final_famd.fit_transform(famd_input)

        # Extract underlying values if it's a DataFrame to avoid index misalignment
        if isinstance(famd_coords, DataFrame):
            famd_values = famd_coords.values
        else:
            famd_values = famd_coords

        # Create DataFrame with proper index and column names
        famd_df = DataFrame(
            famd_values,
            index=famd_input.index,
            columns=[f"FAMD_{i}" for i in range(n_components)],
        ).round(4)

        if debug:
            peek_dataframe(famd_df, "Generated FAMD Components", debug)

        logger.info(
            tabulate(final_famd.eigenvalues_summary, headers="keys", tablefmt="psql")
        )
        logger.info(
            tabulate(final_famd.column_coordinates_, headers="keys", tablefmt="psql")
        )

        return famd_df

    except Exception as e:
        logger.warning(f"FAMD generation failed: {e}")
        return DataFrame(index=famd_input.index)


def main():
    args = parse_args()
    debug = args.debug

    # Configure logging to file
    log_filename = f"{args.cell_type}_{args.modality}_prep_pb.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Logging configured. Writing to {log_filename}")

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
    # if any cell-type specific covariates have missingness fill them, zero; ie missing
    covars_df[counts_term] = covars_df[counts_term].fillna(0)
    if modality == "atac":
        covars_df[probs_term] = covars_df[probs_term].fillna(0)

    if debug:
        logger.debug(covars_df.describe())
        logger.debug(covars_df.info())

    quants_df = load_quantified_data(
        quants_dir, args.project, cell_type, modality, debug
    )

    # scale the dataframe features
    scaled_df = scale_dataframe(quants_df)

    scaled_file = quants_dir / f"{args.project}.{cell_type}.{modality}.scaled.parquet"
    logger.info(f"Saving scaled data to {scaled_file}")
    scaled_df.to_parquet(scaled_file)

    # combine the quantifications and covariates for modeling and cleaning of non-target variance
    data_df = covars_df.merge(quants_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(data_df, "merged the covariates and quantifications", debug)

    # Identify known covariates to save and use for regression before PCA
    known_covariates = [
        "age",
        "bmi",
        "pmi",
        "ph",
        "sex",
        "pool",
        "ancestry",
        "smoker",
        counts_term,
    ]
    if modality == "atac":
        known_covariates.append(probs_term)

    # Generate latent features representing non-target variance base on high variance features
    pca_df = generate_latent_features(
        quants_df,
        data_df,
        known_covariates,
        args.project,
        quants_dir,
        out_figure_path,
        title_suffix,
        args.top_var_fraction,
        args.imputer_type,
        debug,
    )

    # extend the fixed effect terms to include the PCAs
    ext_data_df = data_df.merge(pca_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(ext_data_df, "Extended Data DataFrame", debug)

    # Save final covariates terms to a file
    final_covariates = known_covariates + pca_df.columns.tolist()

    # Generate FAMD features from these final covariates (compressing them)
    famd_df = generate_famd_features(
        ext_data_df,
        [x for x in final_covariates if x != "age"],
        out_figure_path,
        title_suffix,
        debug=debug,
    )

    # Append FAMD features to the extended dataframe
    ext_data_df = ext_data_df.merge(
        famd_df, how="inner", left_index=True, right_index=True
    )

    # Update final covariates list to include FAMD components
    # (keeping original covariates + FAMD components as requested,
    # or should I replace? "appended to them" implies keeping both)
    final_covariates_extended = final_covariates + famd_df.columns.tolist()

    final_covariates_file = (
        info_dir / f"{args.project}.{cell_type}.{modality}.final_covariates.csv"
    )
    logger.info(
        f"Saving final covariates terms (including FAMD) to {final_covariates_file}"
    )

    # Rename columns for standardizing terms in the output file
    rename_map = {counts_term: "cell_counts"}
    if probs_term:
        rename_map[probs_term] = "label_probs"

    ext_data_df[final_covariates_extended].rename(columns=rename_map).to_csv(
        final_covariates_file
    )

    # check if any of the known or generated covariates are correlated with age
    check_correlations(
        ext_data_df[final_covariates_extended],
        "age",
        [x for x in final_covariates_extended if x != "age"],
    )

    # Regress out non-target covariates (everything except age)
    # Using the original uncompressed covariates for regression as per typical workflow
    # unless instructed otherwise. The FAMD is likely for downstream analysis.
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
