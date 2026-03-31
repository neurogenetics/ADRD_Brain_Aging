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
    parser = argparse.ArgumentParser(description="Run mediation analysis (age -> exog -> endo).")
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--endo-modality", type=str, default="rna")
    parser.add_argument("--exog-modality", type=str, default="atac")
    parser.add_argument(
        "--covariates",
        type=str,
        default="all",
        choices=["none", "all", "knowns", "unknowns", "specified"],
    )
    parser.add_argument("--covariates-list", type=str, nargs="+")
    parser.add_argument("--wls-weight-term", type=str, default="cell_counts_endo")
    parser.add_argument("--n-rep", type=int, default=1000, help="Number of bootstrap replications for mediation.")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def run_single_mediation(
    df, endo_feature, exog_feature, exposure, covariates, weight_term, n_rep
):
    try:
        # Construct formulas
        if covariates:
            covars_str = " + ".join(covariates)
            mediator_formula = f'Q("{exog_feature}") ~ {exposure} + {covars_str}'
            outcome_formula = f'Q("{endo_feature}") ~ {exposure} + Q("{exog_feature}") + {covars_str}'
        else:
            mediator_formula = f'Q("{exog_feature}") ~ {exposure}'
            outcome_formula = f'Q("{endo_feature}") ~ {exposure} + Q("{exog_feature}")'

        # Fit models
        mediator_model = smf.wls(mediator_formula, data=df, weights=df[weight_term])
        outcome_model = smf.wls(outcome_formula, data=df, weights=df[weight_term])

        # Mediation analysis
        med = Mediation(
            outcome_model, 
            mediator_model, 
            exposure=exposure, 
            mediator=f'Q("{exog_feature}")'
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
            "prop_mediated_estimate": summary.loc["Prop. mediated (average)", "Estimate"],
            "prop_mediated_pval": summary.loc["Prop. mediated (average)", "P-value"],
            "total_effect_estimate": summary.loc["Total effect", "Estimate"],
            "total_effect_pval": summary.loc["Total effect", "P-value"]
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
            "total_effect_pval": np.nan
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
        c for c in common_bio_cols if c in endo_covars.columns and c in exog_covars.columns
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
    if args.covariates == "none":
        formula_covariates = []
    elif args.covariates == "all":
        formula_covariates = [c for c in covars.columns if c != "age"]
    elif args.covariates == "knowns":
        formula_covariates = [x for x in covars.columns if not ("PCA_" in x) and x != "age"]
    elif args.covariates == "unknowns":
        formula_covariates = [x for x in covars.columns if ("PCA_" in x)]
    elif args.covariates == "specified":
        if not args.covariates_list:
            raise ValueError(
                "--covariates-list must be provided when using --covariates specified"
            )
        formula_covariates = [c for c in args.covariates_list if c != "age"]

    if args.wls_weight_term:
        if args.wls_weight_term in formula_covariates:
            formula_covariates.remove(args.wls_weight_term)
        if args.wls_weight_term not in data_df.columns:
            raise ValueError(
                f"[{cell_type}] wls_weight_term '{args.wls_weight_term}' not found in covariates."
            )
        logger.info(
            f"[{cell_type}] wls_weight_term '{args.wls_weight_term}' will be used in the regression"
        )
    else:
        raise ValueError("A WLS weight term must be provided.")

    logger.info(f"[{cell_type}] Covariates included in formula: {formula_covariates}")

    logger.info(f"[{cell_type}] Starting mediation analysis on {len(ct_sig_results)} pairs...")
    results = []
    
    counter = 0
    total_comparisons = len(ct_sig_results)
    
    for idx, row in ct_sig_results.iterrows():
        endo_term = row["endo_feature"]
        exog_term = row["exog_feature"]
        
        if endo_term not in data_df.columns or exog_term not in data_df.columns:
            continue

        res = run_single_mediation(
            data_df,
            endo_term,
            exog_term,
            exposure="age",
            covariates=formula_covariates,
            weight_term=args.wls_weight_term,
            n_rep=args.n_rep
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
    logger.info(f"Found {total_sig} significantly mediated pairs across all cell types (ACME FDR <= 0.05)")

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