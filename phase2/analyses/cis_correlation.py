import argparse
import concurrent.futures
import logging
from pathlib import Path
import sys

import numpy as np
from pandas import DataFrame as PandasDF
from pandas import read_csv

# Import functions from pseudobulk_regression
from pseudobulk_regression import (
    regression_model,
    compute_fdr,
    load_final_covariates,
    load_quants,
)

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(description="Run cis correlation analysis.")
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
        choices=["ols", "rlm", "glm", "glm_tweedie", "wls", "vwrlm"],
    )
    parser.add_argument("--wls-weight-term", type=str, default="cell_counts")
    parser.add_argument(
        "--covariates",
        type=str,
        default="all",
        choices=["none", "all", "knowns", "unknowns", "specified"],
    )
    parser.add_argument("--covariates-list", type=str, nargs="+")
    parser.add_argument("--endo-modality", type=str, default="rna")
    parser.add_argument("--exog-modality", type=str, default="atac")
    parser.add_argument("--max-dist", type=int, default=1_000_000)
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def cis_regression(
    df: PandasDF,
    endo_term: str,
    exog_term: str,
    formula_covariates: list,
    regression_type: str,
    weight_term: str = None,
    verbose: bool = False,
):
    if len(formula_covariates) > 0:
        covar_term_formula = " + ".join(formula_covariates)
        this_formula = f'Q("{endo_term}") ~ Q("{exog_term}") + {covar_term_formula}'
    else:
        this_formula = f'Q("{endo_term}") ~ Q("{exog_term}")'

    try:
        result = regression_model(
            this_formula,
            df,
            model_type=regression_type,
            weight_term=weight_term,
        )
        ret_exog_term = f'Q("{exog_term}")'
        ret_list = [
            endo_term,
            exog_term,
            result.params["Intercept"],
            result.params[ret_exog_term],
            result.bse[ret_exog_term],
            result.tvalues[ret_exog_term],
            result.pvalues[ret_exog_term],
        ]
        if verbose:
            logger.debug(f"Success for {endo_term} and {exog_term}")
        return ret_list
    except Exception as e:
        if verbose:
            logger.debug(f"Caught error for {endo_term} and {exog_term}: {e}")
        return [endo_term, exog_term] + [np.nan] * 5


def process_cell_type(
    cell_type,
    ct_endo_results,
    ct_exog_results,
    features_df,
    args,
    quants_dir,
    info_dir,
):
    logger.info(f"--- Processing cell type: {cell_type} ---")

    valid_endo = [f for f in ct_endo_results["feature"] if f in features_df.index]
    valid_exog = [f for f in ct_exog_results["feature"] if f in features_df.index]

    endo_features = features_df.loc[valid_endo].copy()
    exo_features = features_df.loc[valid_exog].copy()

    logger.info(
        f"[{cell_type}] Found {len(endo_features)} significant endo features and {len(exo_features)} significant exog features."
    )

    # Find cis-proximal
    logger.info(f"[{cell_type}] Finding cis-proximal pairs...")
    endo_cis_proximal = {}
    exo_unique = set()

    # If the user has duplicates in gene index, handle gracefully (take first match)
    if not endo_features.index.is_unique:
        endo_features = endo_features[~endo_features.index.duplicated(keep="first")]
    if not exo_features.index.is_unique:
        exo_features = exo_features[~exo_features.index.duplicated(keep="first")]

    for chrom in endo_features["chr"].unique():
        chrom_endo = endo_features[endo_features["chr"] == chrom]
        chrom_exo = exo_features[exo_features["chr"] == chrom]
        if chrom_exo.empty:
            continue
        for endo_id, endo_row in chrom_endo.iterrows():
            start_boundary = np.maximum(
                endo_row["start"] - args.max_dist, chrom_exo["start"].min()
            )
            end_boundary = np.minimum(
                endo_row["end"] + args.max_dist, chrom_exo["end"].max()
            )

            found_df = chrom_exo[
                (chrom_exo["start"] >= start_boundary)
                & (chrom_exo["end"] <= end_boundary)
            ]
            if not found_df.empty:
                endo_cis_proximal[endo_id] = found_df.index.tolist()
                exo_unique.update(found_df.index.tolist())

    total_comparisons = sum(len(v) for v in endo_cis_proximal.values())
    logger.info(
        f"[{cell_type}] Identified {len(endo_cis_proximal)} endo features with {total_comparisons} total pairs to test."
    )

    if total_comparisons == 0:
        logger.info(f"[{cell_type}] No valid cis pairs found. Skipping.")
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

    common_bio_cols = [
        "age",
        "bmi",
        "pmi",
        "ph",
        "sex",
        "ancestry",
        "smoker",
        "pool",
    ]
    base_cols = [c for c in common_bio_cols if c in covars.columns]

    # Intersect indices
    common_idx = covars.index.intersection(endo_quants.index).intersection(
        exog_quants.index
    )
    logger.info(f"[{cell_type}] Number of common samples: {len(common_idx)}")
    covars = covars.loc[common_idx]
    endo_quants = endo_quants.loc[common_idx]
    exog_quants = exog_quants.loc[common_idx]

    data_df = covars.copy()

    # Avoid merging entire quants if they are large; subset to what we need
    endo_cols_to_keep = [
        c for c in endo_cis_proximal.keys() if c in endo_quants.columns
    ]
    exog_cols_to_keep = [c for c in exo_unique if c in exog_quants.columns]

    data_df = data_df.join(endo_quants[endo_cols_to_keep])
    # Ignore overlap errors by explicitly picking columns not in data_df
    exog_cols_clean = [c for c in exog_cols_to_keep if c not in data_df.columns]
    data_df = data_df.join(exog_quants[exog_cols_clean])

    # Determine covariates
    if args.covariates == "none":
        formula_covariates = []
    elif args.covariates == "all":
        formula_covariates = covars.columns.tolist()
    elif args.covariates == "knowns":
        formula_covariates = [x for x in covars.columns if not ("PCA_" in x)]
    elif args.covariates == "unknowns":
        formula_covariates = base_cols + [x for x in covars.columns if ("PCA_" in x)]
    elif args.covariates == "specified":
        if not args.covariates_list:
            raise ValueError(
                "--covariates-list must be provided when using --covariates specified"
            )
        formula_covariates = args.covariates_list

    if args.regression_type in ["wls", "vwrlm"] and args.wls_weight_term:
        if args.wls_weight_term in formula_covariates:
            formula_covariates.remove(args.wls_weight_term)
        if args.wls_weight_term not in data_df.columns:
            raise ValueError(
                f"[{cell_type}] wls_weight_term '{args.wls_weight_term}' not found in covariates."
            )
        logger.info(
            f"[{cell_type}] wls_weight_term '{args.wls_weight_term}' will be used in the regression"
        )

    logger.info(f"[{cell_type}] Covariates included in formula: {formula_covariates}")

    # Run regressions
    logger.info(f"[{cell_type}] Starting {args.regression_type} regression...")
    results = []

    counter = 0
    for endo_term, exos in endo_cis_proximal.items():
        if endo_term not in data_df.columns:
            continue
        for exog_term in exos:
            if exog_term not in data_df.columns:
                continue

            res = cis_regression(
                data_df,
                endo_term,
                exog_term,
                formula_covariates,
                args.regression_type,
                args.wls_weight_term,
                args.debug,
            )
            res.append(cell_type)
            results.append(res)

            counter += 1
            if counter % 1000 == 0:
                print(
                    f"[{cell_type}] Processed {counter}/{total_comparisons}...",
                    end="\r",
                )

    print(f"\n[{cell_type}] Completed.")
    return results


def main():
    args = parse_args()

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"

    log_filename = f"{logs_dir}/all_celltypes_{args.endo_modality}_{args.exog_modality}_{args.regression_type}_cis_correlation.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    # Load features
    features_file = quants_dir / f"{args.project}.features.csv"
    logger.info(f"Loading features from {features_file}")
    features_df = read_csv(features_file, index_col=0)

    # Load endo results
    endo_results_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}.all_celltypes.{args.regression_type}_fdr_filtered.age.csv"
    )
    logger.info(f"Loading endo results from {endo_results_file}")
    endo_results = read_csv(endo_results_file)
    endo_results = endo_results[endo_results["fdr_bh"] < args.fdr_threshold]

    # Load exog results
    exog_results_file = (
        # results_dir
        # / f"{args.project}.{args.exog_modality}.all_celltypes.{args.regression_type}d.age.csv"
        results_dir
        / f"{args.project}.all_celltypes.{args.exog_modality}.{args.regression_type}.age.csv"
    )
    logger.info(f"Loading exog results from {exog_results_file}")
    exog_results = read_csv(exog_results_file)
    # exog_results = exog_results[exog_results["fdr_bh"] < args.fdr_threshold]

    cell_types = endo_results["tissue"].unique()

    all_results = []

    logger.info(f"Starting parallel processing for {len(cell_types)} cell types...")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for cell_type in cell_types:
            ct_endo_results = endo_results[endo_results["tissue"] == cell_type]
            ct_exog_results = exog_results[exog_results["tissue"] == cell_type]

            futures.append(
                executor.submit(
                    process_cell_type,
                    cell_type,
                    ct_endo_results,
                    ct_exog_results,
                    features_df,
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
        logger.info("No successful regressions across any cell types. Exiting.")
        return

    results_df = PandasDF(
        data=all_results,
        columns=[
            "endo_feature",
            "exog_feature",
            "intercept",
            "coef",
            "stderr",
            "z",
            "p-value",
            "tissue",
        ],
    )

    results_df["bh_fdr"] = compute_fdr(results_df["p-value"].fillna(1))

    total_sig = results_df.loc[results_df["bh_fdr"] <= 0.05].shape[0]
    logger.info(
        f"Found {total_sig} significant cis pairs across all cell types (FDR <= 0.05)"
    )

    out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.cis.csv"
    )
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved all results to {out_file}")

    sig_out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.cis.fdr_filtered.csv"
    )
    sig_results_df = results_df.loc[results_df["bh_fdr"] <= 0.05]
    sig_results_df.to_csv(sig_out_file, index=False)
    logger.info(f"Saved significant results to {sig_out_file}")

    # Generate and log summary table
    summary_rows = []
    for cell_type in sorted(cell_types):
        ct_endo = endo_results[endo_results["tissue"] == cell_type]["feature"].nunique()
        ct_exog = exog_results[exog_results["tissue"] == cell_type]["feature"].nunique()

        ct_sig = sig_results_df[sig_results_df["tissue"] == cell_type]
        corr_endo = ct_sig["endo_feature"].nunique()
        corr_exog = ct_sig["exog_feature"].nunique()
        pairs = len(ct_sig)

        pct_endo = (corr_endo / ct_endo * 100) if ct_endo > 0 else 0
        pct_exog = (corr_exog / ct_exog * 100) if ct_exog > 0 else 0

        summary_rows.append(
            {
                "Cell Type": cell_type,
                f"Sig {args.endo_modality.upper()}": ct_endo,
                f"Corr {args.endo_modality.upper()}": corr_endo,
                f"% {args.endo_modality.upper()} Corr": f"{pct_endo:.1f}%",
                f"{args.exog_modality.upper()}": ct_exog,
                f"Corr {args.exog_modality.upper()}": corr_exog,
                f"% {args.exog_modality.upper()} Corr": f"{pct_exog:.1f}%",
                "Sig Pairs": pairs,
            }
        )

    summary_df = PandasDF(summary_rows)
    logger.info(
        f"Cis-Correlation Summary (FDR <= 0.05):\n{summary_df.to_string(index=False)}"
    )


if __name__ == "__main__":
    main()
