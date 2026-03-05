import argparse
import logging
from pathlib import Path
import sys

import numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from pandas import DataFrame as PandasDF
from pandas import read_csv, read_parquet
import statsmodels.stats.multitest as smm

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(description="Run pseudobulk regression analysis.")
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
        default="ols",
        choices=["ols", "rlm", "glm", "glm_tweedie", "wls"],
        help="Regression method to use.",
    )
    parser.add_argument(
        "--wls-weight-term",
        type=str,
        default="cell_counts",
        help="Covariate term to be used for weights in WLS regression.",
    )
    parser.add_argument(
        "--covariates",
        type=str,
        default="all",
        choices=["none", "all", "knowns", "unknowns", "specified"],
        help="Which covariates to include. 'none' includes only age. 'all' includes all. 'knowns' excludes PCA_ terms. 'unknowns' includes only age and PCA_ terms. 'specified' requires --covariates-list.",
    )
    parser.add_argument(
        "--covariates-list",
        type=str,
        nargs="+",
        help="List of specific covariate names to use when --covariates specified is selected. 'age' is always included.",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def load_final_covariates(
    info_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> PandasDF:
    final_covars_file = (
        info_dir / f"{project}.{cell_type}.{modality}.final_covariates.csv"
    )
    if not final_covars_file.exists():
        logger.error(f"Final covariates file not found: {final_covars_file}")
        raise FileNotFoundError(f"Final covariates file not found: {final_covars_file}")

    final_covars_df = read_csv(final_covars_file, index_col=0)
    if debug:
        logger.info(f"Loaded final covariates file: {final_covars_file}")
        print(final_covars_df.head())
    return final_covars_df


def load_quants(
    quants_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> PandasDF:
    data_file = quants_dir / f"{project}.{cell_type}.{modality}.parquet"
    if not data_file.exists():
        logger.error(f"quants data file not found: {data_file}")
        raise FileNotFoundError(f"quants data file not found: {data_file}")

    quants_df = read_parquet(data_file)
    if debug:
        logger.info(f"Loaded quants data file: {data_file}")
        print(quants_df.head())
    return quants_df


def regression_model(
    formula: str, df: PandasDF, model_type: str = "ols", weight_term: str = None
):
    if model_type == "ols":
        model = smf.ols(formula=formula, data=df)
        result = model.fit(cov_type="HC3")
    elif model_type == "rlm":
        model = smf.rlm(formula=formula, data=df, M=sm.robust.norms.HuberT())
        result = model.fit()
    elif model_type == "glm_tweedie":
        model = smf.glm(
            formula=formula,
            data=df,
            family=sm.families.Tweedie(
                link=sm.families.links.Log(), var_power=1.6, eql=True
            ),
        )
        result = model.fit()
    elif model_type == "glm":
        model = smf.glm(formula=formula, data=df)
        result = model.fit()
    elif model_type == "wls":
        model = smf.wls(formula=formula, data=df, weights=df[weight_term])
        result = model.fit()
    return result


def age_regression(
    df: PandasDF,
    feature: str,
    regression_type: str,
    verbose: bool = False,
    weight_term: str = None,
    formula_covariates: list = None,
) -> list:
    endo_term = feature
    exog_term = "age"

    if formula_covariates is None:
        covariates = [x for x in df.columns.tolist() if x not in [endo_term, exog_term]]
    else:
        covariates = formula_covariates

    if len(covariates) > 0:
        covar_term_formula = " + ".join(covariates)
        this_formula = f'Q("{endo_term}") ~ {exog_term} + {covar_term_formula}'
    else:
        this_formula = f'Q("{endo_term}") ~ {exog_term}'

    try:
        # run regression via statsmodel
        result = regression_model(
            this_formula,
            df,
            model_type=regression_type,
            weight_term=weight_term,
        )
        ret_list = [
            endo_term,
            result.params["Intercept"],
            result.params[exog_term],
            result.bse[exog_term],
            result.tvalues[exog_term],
            result.pvalues[exog_term],
            np.exp(result.params[exog_term]),
            result.params[exog_term] / np.log(2),
            (np.exp(result.params[exog_term]) - 1) * 100,
        ]
        if verbose:
            print(f"df shape {df.shape}")
            print(result.summary())
            print(
                [
                    "feature",
                    "intercept",
                    "coef",
                    "stderr",
                    "test_statistic",
                    "p-value",
                    "fc",
                    "log2fc",
                    "percentchange",
                ]
            )
            print(ret_list)
    except Exception as e:
        if verbose:
            print(f"Caught Error for {endo_term}: {e}")
        ret_list = [endo_term] + [np.nan] * 5

    return ret_list


def regress_age(
    quants: PandasDF,
    covars: PandasDF,
    cell_name: str,
    modality: str,
    regression_type: str,
    results_dir: Path,
    project: str,
    covariate_type: str,
    covariates_list: list = None,
    wls_weight_term: str = None,
) -> None:
    logger.info(f"Starting regression for {cell_name} using {regression_type}")

    # Merge quants with covariates (we need 'age')
    # Note: covars contains 'age' and other technical factors.
    if covariate_type == "none":
        terms = ["age"]
    elif covariate_type == "all":
        terms = covars.columns.tolist()
    elif covariate_type == "knowns":
        terms = [x for x in covars.columns.tolist() if not x.startswith("PCA_")]
    elif covariate_type == "unknowns":
        terms = ["age"] + [x for x in covars.columns.tolist() if x.startswith("PCA_")]
    elif covariate_type == "specified":
        if not covariates_list:
            raise ValueError(
                "--covariates-list must be provided when using --covariates specified"
            )
        terms = ["age"] + [x for x in covariates_list if x != "age"]
    else:
        raise ValueError(f"Unknown covariate_type: {covariate_type}")

    formula_covariates = [x for x in terms if x != "age"]

    if regression_type == "wls" and wls_weight_term:
        if wls_weight_term in formula_covariates:
            formula_covariates.remove(wls_weight_term)
        if wls_weight_term not in terms:
            terms.append(wls_weight_term)

    # Ensure all terms exist in the dataframe
    missing_terms = [x for x in terms if x not in covars.columns]
    if missing_terms:
        raise ValueError(
            f"The following specified covariates are missing from the data: {missing_terms}"
        )

    logger.info(f"Independent variable and covariates terms used: {terms}")
    data_df = quants.merge(
        covars[terms],
        how="inner",
        left_index=True,
        right_index=True,
    )

    features = quants.columns.tolist()
    logger.info(f"Processing {len(features)} features...")

    type_results = []
    for i, feature in enumerate(features):
        type_results.append(
            age_regression(
                data_df[[feature] + terms],
                feature,
                regression_type,
                weight_term=wls_weight_term,
                formula_covariates=formula_covariates,
            )
        )
        if (i + 1) % 1000 == 0:
            print(f"Processed {i + 1}/{len(features)}...", end="\r")

    print("\nProcessing complete.")

    results_df = PandasDF(
        data=type_results,
        columns=[
            "feature",
            "intercept",
            "coef",
            "stderr",
            "z",
            "p-value",
            "fc",
            "log2fc",
            "percentchange",
        ],
    )
    results_df["bh_fdr_tissue"] = compute_fdr(results_df["p-value"].fillna(1))
    results_df["tissue"] = cell_name

    total_sig = results_df.loc[results_df.bh_fdr_tissue <= 0.05].shape[0]
    logger.info(
        f"in {cell_name} found {total_sig} features that are significant using within tissue FDR"
    )

    # Save results
    out_file = (
        results_dir / f"{project}.{modality}.{cell_name}.{regression_type}.age.csv"
    )
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved results to {out_file}")


def compute_fdr(pvalues):
    bh_adj = smm.fdrcorrection(pvalues)
    return bh_adj[1]


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"

    # Configure logging to file and stdout
    log_filename = f"{logs_dir}/{args.cell_type}_{args.modality}_{args.regression_type}_pseudobulk_regression.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Logging configured. Writing to {log_filename}")

    # Ensure results directory exists
    results_dir.mkdir(parents=True, exist_ok=True)

    modality = args.modality
    cell_type = args.cell_type
    project = args.project

    # Load final covariates (contains 'age' and we might need it,
    # but strictly speaking we just need 'age' which is in the original covariates too.
    # Using final_covariates is safer to ensure sample alignment if any filtering happened).
    covars_df = load_final_covariates(info_dir, project, cell_type, modality, debug)

    # Load quants (clean data)
    quants_df = load_quants(quants_dir, project, cell_type, modality, debug)

    # Run regression
    regress_age(
        quants_df,
        covars_df,
        cell_type,
        modality,
        args.regression_type,
        results_dir,
        project,
        args.covariates,
        args.covariates_list,
        args.wls_weight_term,
    )


if __name__ == "__main__":
    main()
