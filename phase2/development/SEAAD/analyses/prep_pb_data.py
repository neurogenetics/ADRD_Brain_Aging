#!/usr/bin/env python3
"""
Prepare pseudobulk converted data for disease-status regression analysis.
Identifies and regresses out unobserved latent technical/biological noise (PCA)
from autosomal features, while preserving target-disease variance.
"""

import os
import sys
import logging
import argparse
import warnings
from pathlib import Path
import numpy as np
import pandas as pd
from tabulate import tabulate

# Add legacy analyses folder to path to reuse variance_utils helpers
sys.path.append(str(Path(__file__).resolve().parent.parent.parent.parent / "analyses"))

from sklearn.ensemble import ExtraTreesRegressor
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer, KNNImputer, SimpleImputer
from sklearn.linear_model import LinearRegression
from patsy import dmatrix

from variance_utils import (
    get_high_variance_features,
    iterate_model_component_counts,
    component_from_max_curve,
    generate_selected_model,
    check_correlations,
)

# Configure logging
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare disease pseudobulk data and model unobserved non-target variance."
    )
    parser.add_argument(
        "--project",
        type=str,
        default="seaad_ec_multiome",
        help="Project name used for file prefixes.",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad",
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
        required=True,
        help="Cell type to analyze.",
    )
    parser.add_argument(
        "--target-variable",
        type=str,
        default="dx",
        help="The primary disease target variable column in covariates (default: 'dx').",
    )
    parser.add_argument(
        "--known-covariates",
        type=str,
        nargs="+",
        default=["sex", "race", "ageDeath", "PMI", "pH", "brainWeight"],
        help="List of known technical/biological sample-level covariates.",
    )
    parser.add_argument(
        "--imputer-type",
        type=str,
        default="zero",
        choices=["iterative", "knn", "simple", "zero"],
        help="Type of imputation to use.",
    )
    parser.add_argument(
        "--top-var-fraction",
        type=float,
        default=0.10,
        help="Fraction of top variable features to analyze (default: 0.10).",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def peek_dataframe(df: pd.DataFrame, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    logger.info("DataFrame shape: %s", df.shape)
    if verbose:
        if len(df.columns) < 25:
            print(tabulate(df.head(), headers="keys", tablefmt="psql"))
        else:
            print(f"Index head: {df.index.values[0:10]}")
            print(f"Columns head: {df.columns.values[0:10]}")


def load_covariates(
    info_dir: Path, project: str, modality: str, debug: bool = False
) -> pd.DataFrame:
    covars_file = info_dir / f"{project}.covariates.{modality}.csv"
    covars_df = pd.read_csv(covars_file, index_col=0)
    peek_dataframe(covars_df, f"Loaded covariates file: {covars_file}", debug)
    return covars_df


def load_quantified_data(
    quants_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> pd.DataFrame:
    data_file = quants_dir / f"{project}.{cell_type}.{modality}.parquet"
    quants_df = pd.read_parquet(data_file)
    peek_dataframe(quants_df, f"Loaded quantified data file: {data_file}", debug)
    return quants_df


def load_autosomal_features(features_file: Path, debug: bool = False) -> list[str]:
    logger.info("Loading features from %s", features_file)
    features_df = pd.read_csv(features_file)
    autosomes = [f"chr{i}" for i in range(1, 23)]
    autosomal_df = features_df[features_df["chr"].isin(autosomes)]
    if debug:
        logger.debug("Autosomal features shape: %s", autosomal_df.shape)
    return autosomal_df["gene"].tolist()


def perform_regression_correction(
    feature_df: pd.DataFrame,
    covariate_df: pd.DataFrame,
    covariate_cols: list[str],
    debug: bool = False,
) -> pd.DataFrame:
    """
    Regresses out specified covariates from features.
    Handles categorical covariates correctly via statsmodels formula instead of blind one-hot encoding.
    Preserves original NaNs in output residuals.
    """
    logger.info("Regressing out %s from features...", covariate_cols)

    # Align indices
    common_idx = feature_df.index.intersection(covariate_df.index)
    if len(common_idx) < len(feature_df):
        logger.warning(
            "Regression: Dropping %d samples not in covariates.",
            len(feature_df) - len(common_idx),
        )

    Y_orig = feature_df.loc[common_idx]
    X_source = covariate_df.loc[common_idx, covariate_cols].copy()

    if debug:
        peek_dataframe(X_source, "Covariates matrix for regression")

    # Handle missing in Y for fit (impute with mean)
    Y_fit = Y_orig.copy()
    if Y_fit.isnull().values.any():
        logger.warning("Found missing values in features. Imputing with mean for fit.")
        Y_fit = Y_fit.fillna(Y_fit.mean())

    # Build design matrix with Patsy
    formula = "~ " + " + ".join(covariate_cols)
    X_design = dmatrix(formula, X_source, return_type="dataframe")

    # Fit multi-target OLS
    reg = LinearRegression(n_jobs=-1)
    reg.fit(X_design, Y_fit)

    # Residuals = Original - Predicted
    residuals = Y_orig - reg.predict(X_design)
    return residuals


def determine_pca_components(
    imputed_df: pd.DataFrame,
    max_count: int,
    out_prefix: str = None,
    debug: bool = False,
    title_suffix: str = "",
) -> pd.DataFrame:
    logger.info("Determine the number of PCA components to use")

    if max_count <= 1:
        logger.info("Max component count is <= 1. Defaulting to 1 component.")
        num_comp = 1
        pca_mdl, pca_df, _, _ = generate_selected_model(num_comp, imputed_df, "PCA")
        logger.info("PCA DataFrame shape: %s", pca_df.shape)
        peek_dataframe(pca_df, "PCA variance components generated", debug)
        logger.info("Explained variance ratio: %s", pca_mdl.explained_variance_ratio_)
        return pca_df

    r2_values, rmse_values = iterate_model_component_counts(
        max_count, imputed_df, "PCA"
    )
    if debug:
        logger.debug("r2_values: %s", r2_values)
        logger.debug("rmse_values: %s", rmse_values)

    knee_rmse = component_from_max_curve(rmse_values, "RMSE", out_prefix, title_suffix)
    knee_r2 = component_from_max_curve(r2_values, "R2", out_prefix, title_suffix)
    num_comp = max(knee_rmse, knee_r2)
    logger.info("N = %d components will be used", num_comp)

    pca_mdl, pca_df, _, _ = generate_selected_model(num_comp, imputed_df, "PCA")
    logger.info("PCA DataFrame shape: %s", pca_df.shape)
    peek_dataframe(pca_df, "PCA variance components generated", debug)
    logger.info("Explained variance ratio: %s", pca_mdl.explained_variance_ratio_)
    return pca_df


def generate_latent_features(
    quants_df: pd.DataFrame,
    covariates_df: pd.DataFrame,
    covariate_cols: list[str],
    project: str,
    quants_dir: Path,
    out_figure_path: Path,
    title_suffix: str,
    top_var_fraction: float,
    imputer_type: str,
    debug: bool,
) -> pd.DataFrame:
    logger.info("Begin modeling non-target variance in the data")

    # Filter autosomal features
    features_file = quants_dir / f"{project}.features.csv"
    if features_file.exists():
        autosomal_genes = load_autosomal_features(features_file, debug)
        candidate_features = quants_df.columns.intersection(autosomal_genes).tolist()
        logger.info(
            "Restricted to %d autosomal features present in data",
            len(candidate_features),
        )
    else:
        logger.warning(
            "Features file not found at %s. Using all features.", features_file
        )
        candidate_features = quants_df.columns.tolist()

    # Subset quants
    candidate_quants = quants_df[candidate_features]

    # Drop samples with missing values instead of imputing
    clean_quants = candidate_quants.dropna()
    dropped_count = candidate_quants.shape[0] - clean_quants.shape[0]
    logger.info(
        "Dropped %d samples with missing values. Cleaned DataFrame shape: %s",
        dropped_count,
        clean_quants.shape,
    )

    # Regress out known covariates effects before PCA to focus on unknown variance
    residual_df = perform_regression_correction(
        clean_quants, covariates_df, covariate_cols, debug
    )

    # High variance selection on residuals
    variance_features = get_high_variance_features(residual_df, top_var_fraction)
    logger.info("Found %d high variance features in residuals", len(variance_features))
    max_count = int(
        min(
            residual_df[variance_features].shape[0],
            residual_df[variance_features].shape[1],
        )
        / 2
    )
    logger.info("Max components count: %d", max_count)

    # Perform PCA on residuals
    pca_df = determine_pca_components(
        residual_df[variance_features],
        max_count,
        str(out_figure_path),
        debug,
        title_suffix,
    )
    return pca_df


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    figs_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"

    # Configure logging to file
    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    log_filename = logs_dir / f"{safe_ct}_{args.modality}_disease_prep_pb.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
        force=True,
    )
    logger.info("Command line: %s", " ".join(sys.argv))
    logger.info("Logging configured. Writing to %s", log_filename)

    figs_dir.mkdir(parents=True, exist_ok=True)

    modality = args.modality
    cell_type = args.cell_type
    counts_term = f"{safe_ct}_counts"

    out_figure_path = figs_dir / f"{args.project}_{safe_ct}_{modality}"
    title_suffix = f"{cell_type} ({modality.upper()})"

    # Load covariates
    covars_df = load_covariates(info_dir, args.project, modality, debug)

    # Filter out samples where the target variable is missing/NaN (not part of study)
    if args.target_variable not in covars_df.columns:
        logger.error(
            "Target variable '%s' not found in loaded covariates.", args.target_variable
        )
        sys.exit(1)

    initial_len = len(covars_df)
    covars_df = covars_df.dropna(subset=[args.target_variable])
    filtered_len = len(covars_df)
    if initial_len != filtered_len:
        logger.info(
            "Dropped %d samples missing target variable '%s'. Remaining: %d",
            initial_len - filtered_len,
            args.target_variable,
            filtered_len,
        )

    # Ensure count term is filled
    covars_df[counts_term] = covars_df[counts_term].fillna(0)

    # Load quantified parquet files
    quants_df = load_quantified_data(
        quants_dir, args.project, cell_type, modality, debug
    )

    # Merge covariates and quantifications
    data_df = covars_df.merge(quants_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(data_df, "merged covariates and quantifications", debug)

    # Build known covariates list
    known_covariates = [args.target_variable] + args.known_covariates + [counts_term]

    # Model unobserved non-target variance regressing out the primary target variable first
    pca_df = generate_latent_features(
        quants_df,
        data_df,
        [args.target_variable],
        args.project,
        quants_dir,
        out_figure_path,
        title_suffix,
        args.top_var_fraction,
        args.imputer_type,
        debug,
    )

    # Merge covariates with PCA components
    ext_data_df = data_df.merge(pca_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(ext_data_df, "Extended Data DataFrame", debug)

    # Prepare list for final covariates output
    final_covariates = known_covariates + pca_df.columns.tolist()

    # Save final covariates file
    final_covariates_file = (
        info_dir / f"{args.project}.{safe_ct}.{modality}.final_covariates.csv"
    )
    logger.info("Saving final covariates terms to %s", final_covariates_file)

    # Rename cell type counts column to standardized term 'cell_counts'
    rename_map = {counts_term: "cell_counts"}
    ext_data_df[final_covariates].rename(columns=rename_map).to_csv(
        final_covariates_file
    )

    # Perform covariate correlations check against the target variable
    logger.info(
        "Checking for correlations between %s and known covariates",
        args.target_variable,
    )
    check_correlations(
        ext_data_df[known_covariates],
        args.target_variable,
        [x for x in known_covariates if x != args.target_variable],
    )
    logger.info(
        "Checking for correlations between %s and PCA covariates", args.target_variable
    )
    check_correlations(
        ext_data_df[final_covariates],
        args.target_variable,
        [x for x in final_covariates if x.startswith("PCA_")],
    )

    # Regress out PCA unobserved factors from quantifications to compute residual parquet
    feature_cols = quants_df.columns.tolist()
    residuals_df = perform_regression_correction(
        ext_data_df[feature_cols],
        ext_data_df,
        pca_df.columns.tolist(),
        debug,
    )

    # Save clean residuals parquet
    residuals_file = (
        quants_dir / f"{args.project}.{safe_ct}.{modality}.residuals.parquet"
    )
    logger.info("Saving clean residuals to %s", residuals_file)
    residuals_df.to_parquet(residuals_file)

    logger.info("=== Disease Pseudobulk Prep Finished Successfully ===")


if __name__ == "__main__":
    main()
