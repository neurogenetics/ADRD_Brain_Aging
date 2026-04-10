import argparse
import concurrent.futures
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.mediation import Mediation

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
        description="Run mediation analysis (age -> exog -> endo)."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--endo-modality", type=str, default="rna")
    parser.add_argument("--exog-modality", type=str, default="atac")

    parser.add_argument(
        "--endo-covariates",
        type=str,
        default="all",
        choices=["none", "all", "knowns", "unknowns", "specified"],
    )
    parser.add_argument("--endo-covariates-list", type=str, nargs="+")
    parser.add_argument("--endo-weight-term", type=str, default="cell_counts_endo")

    parser.add_argument(
        "--exog-covariates",
        type=str,
        default="all",
        choices=["none", "all", "knowns", "unknowns", "specified"],
    )
    parser.add_argument("--exog-covariates-list", type=str, nargs="+")
    parser.add_argument("--exog-weight-term", type=str, default="cell_counts_exog")

    parser.add_argument(
        "--n-rep",
        type=int,
        default=1000,
        help="Number of bootstrap replications for mediation.",
    )
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def run_single_mediation(
    df,
    endo_feature,
    exog_feature,
    exposure,
    endo_covariates,
    exog_covariates,
    endo_weight_term,
    exog_weight_term,
    n_rep,
):
    try:
        # Construct formulas
        exog_covars_str = " + ".join(exog_covariates) if exog_covariates else ""
        mediator_formula = f'Q("{exog_feature}") ~ {exposure}'
        if exog_covars_str:
            mediator_formula += f" + {exog_covars_str}"

        endog_covars_str = " + ".join(endo_covariates) if endo_covariates else ""
        outcome_formula = f'Q("{endo_feature}") ~ {exposure} + Q("{exog_feature}")'
        if endog_covars_str:
            outcome_formula += f" + {endog_covars_str}"

        # Fit models
        mediator_model = smf.wls(
            mediator_formula, data=df, weights=df[exog_weight_term]
        )
        outcome_model = smf.wls(outcome_formula, data=df, weights=df[endo_weight_term])

        # Mediation analysis
        med = Mediation(
            outcome_model,
            mediator_model,
            exposure=exposure,
            mediator=f'Q("{exog_feature}")',
        )
        med_result = med.fit(n_rep=n_rep)
        summary = med_result.summary()

        # We extract ACME (average), ADE (average), Prop. mediated (average)
        return {
            "endo_feature": endo_feature,
            "exog_feature": exog_feature,
            "acme_estimate": summary.loc["ACME (average)", "Estimate"],
            "acme_pval": summary.loc["ACME (average)", "P-value"],
            "ade_estimate": summary.loc["ADE (average)", "Estimate"],
            "ade_pval": summary.loc["ADE (average)", "P-value"],
            "prop_mediated_estimate": summary.loc[
                "Prop. mediated (average)", "Estimate"
            ],
            "prop_mediated_pval": summary.loc["Prop. mediated (average)", "P-value"],
            "total_effect_estimate": summary.loc["Total effect", "Estimate"],
            "total_effect_pval": summary.loc["Total effect", "P-value"],
        }
    except Exception as e:
        logger.debug(f"Mediation failed for {endo_feature} - {exog_feature}: {e}")
        return {
            "endo_feature": endo_feature,
            "exog_feature": exog_feature,
            "acme_estimate": np.nan,
            "acme_pval": np.nan,
            "ade_estimate": np.nan,
            "ade_pval": np.nan,
            "prop_mediated_estimate": np.nan,
            "prop_mediated_pval": np.nan,
            "total_effect_estimate": np.nan,
            "total_effect_pval": np.nan,
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
    exog_covars = load_final_covariates(
        info_dir, args.project, cell_type, args.exog_modality, args.debug
    )

    # Merge covariates
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
    base_cols = [
        c
        for c in common_bio_cols
        if c in endo_covars.columns and c in exog_covars.columns
    ]

    endo_spec = [c for c in endo_covars.columns if c not in base_cols]
    exog_spec = [c for c in exog_covars.columns if c not in base_cols]

    covars = endo_covars[base_cols].copy()
    for c in endo_spec:
        covars[f"{c}_endo"] = endo_covars[c]
    for c in exog_spec:
        covars[f"{c}_exog"] = exog_covars[c]

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
    def resolve_covariates(arg_covar, arg_covar_list, original_cols, spec_cols, suffix):
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

        resolved = []
        for c in raw_covars:
            if c in base_cols:
                resolved.append(c)
            elif c in spec_cols:
                resolved.append(f"{c}_{suffix}")
            else:
                resolved.append(c)
        return list(dict.fromkeys(resolved))

    endo_formula_covariates = resolve_covariates(
        args.endo_covariates,
        args.endo_covariates_list,
        endo_covars.columns,
        endo_spec,
        "endo",
    )
    exog_formula_covariates = resolve_covariates(
        args.exog_covariates,
        args.exog_covariates_list,
        exog_covars.columns,
        exog_spec,
        "exog",
    )

    endo_weight_term = args.endo_weight_term
    exog_weight_term = args.exog_weight_term

    if endo_weight_term in endo_formula_covariates:
        endo_formula_covariates.remove(endo_weight_term)
    if endo_weight_term in exog_formula_covariates:
        exog_formula_covariates.remove(endo_weight_term)

    if exog_weight_term in exog_formula_covariates:
        exog_formula_covariates.remove(exog_weight_term)
    if exog_weight_term in endo_formula_covariates:
        endo_formula_covariates.remove(exog_weight_term)

    if endo_weight_term not in data_df.columns:
        raise ValueError(
            f"[{cell_type}] endo_weight_term '{endo_weight_term}' not found in covariates."
        )
    if exog_weight_term not in data_df.columns:
        raise ValueError(
            f"[{cell_type}] exog_weight_term '{exog_weight_term}' not found in covariates."
        )

    all_covar_cols = list(
        dict.fromkeys(
            ["age", endo_weight_term, exog_weight_term]
            + endo_formula_covariates
            + exog_formula_covariates
        )
    )

    # Check consistency globally across all expected covariates to ensure identically matched samples
    # before evaluating individual genes. Features with missing values aren't an issue unless they
    # are specifically evaluated in a pair, which we'll handle gracefully.
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
        f"[{cell_type}] Exog covariates: {exog_formula_covariates} | Weight: {exog_weight_term}"
    )

    logger.info(
        f"[{cell_type}] Starting mediation analysis on {len(ct_sig_results)} pairs..."
    )
    results = []

    counter = 0
    total_comparisons = len(ct_sig_results)

    for idx, row in ct_sig_results.iterrows():
        endo_term = row["endo_feature"]
        exog_term = row["exog_feature"]

        if endo_term not in data_df.columns or exog_term not in data_df.columns:
            continue

        # Isolate exactly the columns we need for this pair and drop missing sample values
        # just for these specific quantitative features (already pre-cleared for covariates above)
        pair_cols = all_covar_cols + [endo_term, exog_term]
        model_df = data_df[pair_cols].dropna()

        if len(model_df) < 10:
            logger.debug(
                f"[{cell_type}] Mediation failed for {endo_term} - {exog_term}: Too few samples ({len(model_df)})."
            )
            continue

        res = run_single_mediation(
            model_df,
            endo_term,
            exog_term,
            exposure="age",
            endo_covariates=endo_formula_covariates,
            exog_covariates=exog_formula_covariates,
            endo_weight_term=endo_weight_term,
            exog_weight_term=exog_weight_term,
            n_rep=args.n_rep,
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

    log_filename = f"{logs_dir}/all_celltypes_{args.endo_modality}_{args.exog_modality}_wls_mediation.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    # Load cis fdr filtered results
    # Defaulting to wls regression type for the input file name since we only do wls mediation
    input_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.wls.cis.fdr_filtered.csv"
    )
    logger.info(f"Loading significant cis correlations from {input_file}")

    if not input_file.exists():
        logger.error(f"Input file not found: {input_file}")
        sys.exit(1)

    sig_results = pd.read_csv(input_file)
    logger.info(f"Loaded {len(sig_results)} significant pairs.")

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
        logger.info("No successful mediations across any cell types. Exiting.")
        return

    results_df = pd.DataFrame(all_results)

    # Filter out missing results
    results_df = results_df.dropna(subset=["acme_estimate"])

    if len(results_df) == 0:
        logger.info("All mediation runs failed or resulted in NaNs. Exiting.")
        return

    # Calculate FDR on ACME (Average Causal Mediation Effect) p-value
    results_df["acme_fdr"] = compute_fdr(results_df["acme_pval"])

    total_sig = results_df.loc[results_df["acme_fdr"] <= 0.05].shape[0]
    logger.info(
        f"Found {total_sig} significantly mediated pairs across all cell types (ACME FDR <= 0.05)"
    )

    out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.wls.mediation.csv"
    )
    results_df.to_csv(out_file, index=False)
    logger.info(f"Saved mediation results to {out_file}")

    sig_out_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.wls.mediation.fdr_filtered.csv"
    )
    sig_results_df = results_df.loc[results_df["acme_fdr"] <= 0.05]
    sig_results_df.to_csv(sig_out_file, index=False)
    logger.info(f"Saved significant mediation results to {sig_out_file}")


if __name__ == "__main__":
    main()
