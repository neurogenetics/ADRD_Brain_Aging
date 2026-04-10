import argparse
import concurrent.futures
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

# Import functions from pseudobulk_regression
from pseudobulk_regression import (
    compute_fdr,
    load_final_covariates,
    load_quants,
)

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run conditioned regression analysis (endo ~ age + exog + covariates)."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--endo-modality", type=str, default="rna")
    parser.add_argument("--exog-modality", type=str, default="atac")

    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
    )

    parser.add_argument(
        "--endo-covariates",
        type=str,
        default="all",
        choices=["none", "all", "knowns", "unknowns", "specified"],
    )
    parser.add_argument("--endo-covariates-list", type=str, nargs="+")
    parser.add_argument("--endo-weight-term", type=str, default="cell_counts")

    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--debug", action="store_true")

    return parser.parse_args()


def run_single_conditioned_regression(
    df,
    endo_feature,
    exog_feature,
    exposure,
    endo_covariates,
    endo_weight_term,
):
    try:
        # Construct formula
        outcome_covars_list = list(dict.fromkeys(endo_covariates))
        outcome_covars_str = (
            " + ".join(outcome_covars_list) if outcome_covars_list else ""
        )
        outcome_formula = f'Q("{endo_feature}") ~ {exposure} + Q("{exog_feature}")'
        if outcome_covars_str:
            outcome_formula += f" + {outcome_covars_str}"

        # Fit model using WLS, similar to the outcome_model in mediation
        outcome_model = smf.wls(outcome_formula, data=df, weights=df[endo_weight_term])
        result = outcome_model.fit()

        return {
            "endo_feature": endo_feature,
            "exog_feature": exog_feature,
            "exposure_coef": result.params[exposure],
            "exposure_stderr": result.bse[exposure],
            "exposure_tval": result.tvalues[exposure],
            "exposure_pval": result.pvalues[exposure],
        }
    except Exception as e:
        logger.debug(
            f"Conditioned regression failed for {endo_feature} - {exog_feature}: {e}"
        )
        return {
            "endo_feature": endo_feature,
            "exog_feature": exog_feature,
            "exposure_coef": np.nan,
            "exposure_stderr": np.nan,
            "exposure_tval": np.nan,
            "exposure_pval": np.nan,
        }


def process_cell_type(
    cell_type,
    ct_sig_results,
    args,
    quants_dir,
    info_dir,
):
    logger.info(f"--- Processing cell type: {cell_type} ---")

    unique_endo = set(ct_sig_results["endo_feature"])
    unique_exog = set(ct_sig_results["exog_feature"])

    if len(ct_sig_results) == 0:
        logger.info(f"[{cell_type}] No significant cis pairs to process.")
        return []

    # Load quants
    logger.info(f"[{cell_type}] Loading quantifications...")
    endo_quants = load_quants(
        quants_dir, args.project, cell_type, args.endo_modality, args.debug
    )
    exog_quants = load_quants(
        quants_dir, args.project, cell_type, args.exog_modality, args.debug
    )

    # Load covariates
    logger.info(f"[{cell_type}] Loading covariates...")
    endo_covars = load_final_covariates(
        info_dir, args.project, cell_type, args.endo_modality, args.debug
    )

    covars = endo_covars.copy()

    # Intersect indices
    common_idx = covars.index.intersection(endo_quants.index).intersection(
        exog_quants.index
    )
    logger.info(f"[{cell_type}] Number of common samples: {len(common_idx)}")
    covars = covars.loc[common_idx]
    endo_quants = endo_quants.loc[common_idx]
    exog_quants = exog_quants.loc[common_idx]

    data_df = covars.copy()

    endo_cols_to_keep = [c for c in unique_endo if c in endo_quants.columns]
    exog_cols_to_keep = [c for c in unique_exog if c in exog_quants.columns]

    data_df = data_df.join(endo_quants[endo_cols_to_keep])
    exog_cols_clean = [c for c in exog_cols_to_keep if c not in data_df.columns]
    data_df = data_df.join(exog_quants[exog_cols_clean])

    # Determine covariates
    def resolve_covariates(arg_covar, arg_covar_list, original_cols):
        if arg_covar == "none":
            raw_covars = []
        elif arg_covar == "all":
            raw_covars = [c for c in original_cols if c != "age"]
        elif arg_covar == "knowns":
            raw_covars = [c for c in original_cols if not ("PCA_" in c) and c != "age"]
        elif arg_covar == "unknowns":
            raw_covars = [c for c in original_cols if ("PCA_" in c)]
        elif arg_covar == "specified":
            raw_covars = [c for c in arg_covar_list if c != "age"]
        else:
            raw_covars = []

        return list(dict.fromkeys(raw_covars))

    endo_formula_covariates = resolve_covariates(
        args.endo_covariates,
        args.endo_covariates_list,
        covars.columns,
    )

    endo_weight_term = args.endo_weight_term

    if endo_weight_term in endo_formula_covariates:
        endo_formula_covariates.remove(endo_weight_term)

    if endo_weight_term not in data_df.columns:
        raise ValueError(
            f"[{cell_type}] endo_weight_term '{endo_weight_term}' not found in covariates."
        )

    all_covar_cols = list(
        dict.fromkeys(["age", endo_weight_term] + endo_formula_covariates)
    )

    before_drop = len(data_df)
    data_df = data_df.dropna(subset=all_covar_cols)
    after_drop = len(data_df)

    logger.info(
        f"[{cell_type}] Dropped {before_drop - after_drop} samples due to missing covariate/weight values. Remaining: {after_drop}"
    )

    logger.info(
        f"[{cell_type}] Endo covariates: {endo_formula_covariates} | Weight: {endo_weight_term}"
    )

    logger.info(
        f"[{cell_type}] Starting conditioned regression analysis on {len(ct_sig_results)} pairs..."
    )
    results = []

    counter = 0
    total_comparisons = len(ct_sig_results)

    for idx, row in ct_sig_results.iterrows():
        endo_term = row["endo_feature"]
        exog_term = row["exog_feature"]

        if endo_term not in data_df.columns or exog_term not in data_df.columns:
            continue

        pair_cols = all_covar_cols + [endo_term, exog_term]
        model_df = data_df[pair_cols].dropna()

        if len(model_df) < 10:
            logger.debug(
                f"[{cell_type}] Conditioned regression failed for {endo_term} - {exog_term}: Too few samples ({len(model_df)})."
            )
            continue

        res = run_single_conditioned_regression(
            model_df,
            endo_term,
            exog_term,
            exposure="age",
            endo_covariates=endo_formula_covariates,
            endo_weight_term=endo_weight_term,
        )
        res["tissue"] = cell_type
        results.append(res)

        counter += 1
        if counter % 100 == 0:
            print(f"[{cell_type}] Processed {counter}/{total_comparisons}...", end="\r")

    print(f"\n[{cell_type}] Completed.")
    return results


def main():
    args = parse_args()

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"

    log_filename = f"{logs_dir}/all_celltypes_{args.endo_modality}_{args.exog_modality}_{args.regression_type}_conditioned_regression.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    # Construct file paths for endo, exog, and cis raw results
    endo_results_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}.all_celltypes.{args.regression_type}_fdr_filtered.age.csv"
    )
    exog_results_file = (
        results_dir
        / f"{args.project}.{args.exog_modality}.all_celltypes.{args.regression_type}_fdr_filtered.age.csv"
    )
    cis_results_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.cis.csv"
    )

    for f in [endo_results_file, exog_results_file, cis_results_file]:
        if not f.exists():
            logger.error(f"Required file not found: {f}")
            sys.exit(1)

    logger.info(f"Loading endo results from {endo_results_file}")
    endo_results = pd.read_csv(endo_results_file)
    endo_results = endo_results[endo_results["fdr_bh"] < args.fdr_threshold]

    logger.info(f"Loading exog results from {exog_results_file}")
    exog_results = pd.read_csv(exog_results_file)
    exog_results = exog_results[exog_results["fdr_bh"] < args.fdr_threshold]

    logger.info(f"Loading cis results from {cis_results_file}")
    results_df = pd.read_csv(cis_results_file)

    # Restrict results to just age associated features in both modalities
    results_df = results_df[
        (results_df["endo_feature"].isin(endo_results["feature"]))
        & (results_df["exog_feature"].isin(exog_results["feature"]))
    ]

    # Recompute the B&H FDR for just the age associated results
    # To use `compute_fdr` from `pseudobulk_regression.py` on the p-values
    results_df["bh_fdr"] = compute_fdr(results_df["p-value"].fillna(1))
    sig_results = results_df[results_df["bh_fdr"] < args.fdr_threshold]

    logger.info(f"Loaded {len(sig_results)} significant age-associated cis pairs.")

    cell_types = sig_results["tissue"].unique()
    all_results = []

    logger.info(f"Starting parallel processing for {len(cell_types)} cell types...")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for cell_type in cell_types:
            ct_sig_results = sig_results[sig_results["tissue"] == cell_type]

            futures.append(
                executor.submit(
                    process_cell_type,
                    cell_type,
                    ct_sig_results,
                    args,
                    quants_dir,
                    info_dir,
                )
            )

        for future in concurrent.futures.as_completed(futures):
            try:
                ct_results = future.result()
                if ct_results:
                    all_results.extend(ct_results)
            except Exception as e:
                logger.error(f"Exception during parallel processing: {e}")

    print("\nAll processing complete.")

    if not all_results:
        logger.info(
            "No successful conditioned regressions across any cell types. Exiting."
        )
        return

    final_results_df = pd.DataFrame(all_results)

    # Filter out missing results
    final_results_df = final_results_df.dropna(subset=["exposure_coef"])

    if len(final_results_df) == 0:
        logger.info("All runs failed or resulted in NaNs. Exiting.")
        return

    # Calculate FDR on the conditioned regression exposure p-value
    final_results_df["exposure_fdr"] = compute_fdr(final_results_df["exposure_pval"])

    total_tested = len(final_results_df)
    nominally_sig = final_results_df.loc[
        final_results_df["exposure_pval"] <= 0.05
    ].shape[0]
    total_sig = final_results_df.loc[
        final_results_df["exposure_fdr"] <= args.fdr_threshold
    ].shape[0]

    logger.info(f"Total pairs tested: {total_tested}")
    logger.info(
        f"Total nominally significant conditioned pairs (p-value <= 0.05): {nominally_sig}"
    )
    logger.info(
        f"Found {total_sig} significantly conditioned pairs across all cell types (Exposure FDR <= {args.fdr_threshold})"
    )

    out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.conditioned.csv"
    )
    final_results_df.to_csv(out_file, index=False)
    logger.info(f"Saved conditioned regression results to {out_file}")

    sig_out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.conditioned.fdr_filtered.csv"
    )
    sig_results_df = final_results_df.loc[
        final_results_df["exposure_fdr"] <= args.fdr_threshold
    ]
    sig_results_df.to_csv(sig_out_file, index=False)
    logger.info(f"Saved significant conditioned regression results to {sig_out_file}")


if __name__ == "__main__":
    main()
