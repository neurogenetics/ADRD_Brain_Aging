import argparse
import logging
from pathlib import Path
from pandas import DataFrame, read_csv, read_parquet
import statsmodels.formula.api as smf
import concurrent.futures
import seaborn as sns
import matplotlib.pyplot as plt
from tabulate import tabulate
from variance_utils import perform_variance_partition, check_correlations, fit_and_report_correlation

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run variance partition analysis on single-cell/pseudobulk data."
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
    if not covars_file.exists():
        logger.error(f"Covariates file not found: {covars_file}")
        raise FileNotFoundError(f"Covariates file not found: {covars_file}")
    
    covars_df = read_csv(covars_file, index_col=0)
    peek_dataframe(covars_df, f"Loaded covariates file: {covars_file}", debug)
    return covars_df


def load_final_covariates(
    info_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> DataFrame:
    final_covars_file = info_dir / f"{project}.{cell_type}.{modality}.final_covariates.csv"
    if not final_covars_file.exists():
        logger.warning(f"Final covariates file not found: {final_covars_file}. skipping final analysis.")
        return None

    final_covars_df = read_csv(final_covars_file, index_col=0)
    peek_dataframe(final_covars_df, f"Loaded final covariates file: {final_covars_file}", debug)
    return final_covars_df


def load_quantified_data(
    quants_dir: Path, project: str, cell_type: str, modality: str, debug: bool = False
) -> DataFrame:
    data_file = quants_dir / f"{project}.{cell_type}.{modality}.parquet"
    if not data_file.exists():
        logger.error(f"Quantified data file not found: {data_file}")
        raise FileNotFoundError(f"Quantified data file not found: {data_file}")

    quants_df = read_parquet(data_file)
    peek_dataframe(quants_df, f"Loaded quantified data file: {data_file}", debug)
    return quants_df


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


def run_variance_partition_parallel(
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


def main():
    args = parse_args()
    debug = args.debug

    # Configure logging to file
    log_filename = f"{args.cell_type}_{args.modality}_variance_partition.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler()
        ],
        force=True
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
    quants_df = load_quantified_data(
        quants_dir, args.project, cell_type, modality, debug
    )
    
    # Merge for Initial Analysis
    data_df = covars_df.merge(quants_df, how="inner", left_index=True, right_index=True)
    peek_dataframe(data_df, "merged the covariates and quantifications", debug)

    # --- Initial Analysis (Known Covariates Only) ---
    logger.info("--- Starting Initial Variance Partition (Known Covariates) ---")
    
    check_covariates = [
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
        check_covariates.append(probs_term)

    check_correlations(data_df, "age", check_covariates)

    fixed_effects = ["age", "bmi", "pmi", "ph", counts_term]
    if modality == "atac":
        fixed_effects.append(probs_term)

    random_effects = ["sex", "pool", "ancestry", "smoker"]

    results = run_variance_partition_parallel(
        data_df, quants_df.columns.tolist(), fixed_effects, random_effects, debug
    )

    process_variance_results(results, out_figure_path, "_known", debug, title_suffix)

    # --- Final Analysis (Known Covariates + PCA) ---
    final_covars_df = load_final_covariates(info_dir, args.project, cell_type, modality, debug)
    
    if final_covars_df is not None:
        logger.info("--- Starting Final Variance Partition (Known + PCA) ---")
        
        # Merge final covariates with quant data
        # Note: final_covars_df likely includes known covariates + PCA + possibly other columns
        # We need to make sure we don't duplicate columns if we merge
        
        # Let's see what final_covars_df has. It should have index as samples.
        # We can just merge final_covars_df with quants_df. 
        # But final_covars_df might already have the known covariates in it?
        # In prep_pb_data.py, final_covariates file was saved from `ext_data_df[final_covariates]`.
        # `ext_data_df` was `data_df` merged with `pca_df`.
        # `data_df` was `covars_df` merged with `quants_df`.
        # So `final_covariates.csv` contains known covariates + PCA columns.
        
        # So if we load it, we have everything except quantifications.
        
        ext_data_df = final_covars_df.merge(quants_df, how="inner", left_index=True, right_index=True)
        peek_dataframe(ext_data_df, "Extended Data DataFrame (Final)", debug)

        # Identify PCA columns
        # They typically start with PCA_ or Factor_ ?
        # In prep_pb_data.py, they are generated by `generate_selected_model` which adds prefix f"{model_type}_".
        # Default is PCA, so "PCA_1", "PCA_2", etc.
        
        pca_cols = [c for c in final_covars_df.columns if c.startswith("PCA_")]
        if not pca_cols:
             # Fallback if NMF or ICA was used, though prep_pb_data defaults to PCA
             pca_cols = [c for c in final_covars_df.columns if c.startswith("NMF_") or c.startswith("ICA_")]
        
        logger.info(f"Identified {len(pca_cols)} PCA components: {pca_cols}")
        
        # Check Age vs PCA correlations
        check_correlations(
            ext_data_df, "age", pca_cols, "check age vs PCA correlations"
        )
        
        # Check PCAs against known covariates
        # We need to know which are "known covariates".
        # We can reuse the lists `fixed_effects` and `random_effects` defined above.
        known_covariates = [x for x in (fixed_effects + random_effects) if x in ext_data_df.columns]
        
        # Note: In prep_pb_data.py, counts_term was renamed to "cell_counts" and probs_term to "label_probs"
        # in the final covariates file.
        # So we need to adjust our known_covariates list to match the column names in final_covars_df.
        
        rename_map_inv = {
            "cell_counts": counts_term,
            "label_probs": probs_term
        }
        # Actually, we should probably check if "cell_counts" exists and use it instead of counts_term
        
        if "cell_counts" in ext_data_df.columns and counts_term not in ext_data_df.columns:
            logger.info(f"Using 'cell_counts' instead of '{counts_term}'")
            if counts_term in known_covariates:
                known_covariates.remove(counts_term)
                known_covariates.append("cell_counts")
            if counts_term in fixed_effects:
                fixed_effects.remove(counts_term)
                fixed_effects.append("cell_counts")

        if probs_term and "label_probs" in ext_data_df.columns and probs_term not in ext_data_df.columns:
            logger.info(f"Using 'label_probs' instead of '{probs_term}'")
            if probs_term in known_covariates:
                known_covariates.remove(probs_term)
                known_covariates.append("label_probs")
            if probs_term in fixed_effects:
                fixed_effects.remove(probs_term)
                fixed_effects.append("label_probs")

        check_pcas_against_known_covariates(
            ext_data_df,
            pca_cols,
            known_covariates,
            out_figure_path,
            title_suffix,
        )
        
        # Prepare lists for variance partition
        final_fixed_effects = fixed_effects + pca_cols
        final_random_effects = random_effects
        
        results_final = run_variance_partition_parallel(
            ext_data_df, quants_df.columns.tolist(), final_fixed_effects, final_random_effects, debug
        )

        process_variance_results(results_final, out_figure_path, "_final", debug, title_suffix)


if __name__ == "__main__":
    main()
