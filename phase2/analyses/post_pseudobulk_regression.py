import argparse
import logging
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas import DataFrame, concat, read_csv
from statsmodels.stats.multitest import multipletests

# Configure logging
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run post-processing for pseudobulk regression analysis."
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
        choices=["rna", "atac"],
        help="Data modality (rna or atac).",
    )
    parser.add_argument(
        "--regression-type",
        type=str,
        default="glm_tweedie",
        choices=["glm", "glm_tweedie", "rlm"],
        help="Regression method to use.",
    )
    parser.add_argument(
        "--min-rlm-effect",
        type=float,
        default=0.001,
        help="Minimum effect size for RLM significant results.",
    )
    parser.add_argument(
        "--debug", 
        action="store_true", 
        help="Enable debug output."
    )
    return parser.parse_args()


def compute_bh_fdr(df: DataFrame, alpha: float = 0.05, p_col: str = 'p-value',
                   method: str = 'fdr_bh', verbose: bool = True) -> DataFrame:
    ret_df = df.copy()
    # Handle NaNs in p-values by filling with 1 before correction
    p_vals = ret_df[p_col].fillna(1)
    test_adjust = multipletests(np.array(p_vals), alpha=alpha, method=method)
    ret_df[method] = test_adjust[1]
    if verbose:
        logger.info(f'Total significant after correction: {ret_df.loc[ret_df[method] < alpha].shape[0]}')
    return ret_df


def volcano_plot(df: DataFrame, project: str, modality: str, regression_type: str, 
                 figures_dir: Path, x_term: str = 'coef', y_term: str = 'p-value', 
                 alpha: float = 0.05, adj_p_col: str = 'fdr_bh', title: str = None, 
                 filter_nseeff: bool = True, extreme_size: float = 10.0):
    
    plot_df = df.copy()
    plot_df = plot_df.reset_index(drop=True)    
    
    # Filter extreme values for plotting
    if filter_nseeff:
        plot_df = plot_df.loc[
            ((-extreme_size < plot_df[x_term]) & 
             (plot_df[x_term] < extreme_size) &
             (~plot_df['z'].isna()) | 
             (plot_df[adj_p_col] < alpha))
        ]
        
    plt.figure(figsize=(9, 9))
    
    # Avoid log(0)
    plot_df[y_term] = plot_df[y_term].replace(0, np.finfo(float).eps)
    log_pvalue = -np.log10(plot_df[y_term])
    
    is_sig = plot_df[adj_p_col] < alpha
    
    sns.set_style("whitegrid")
    sns.scatterplot(x=plot_df[x_term], y=log_pvalue, hue=is_sig, palette={True: 'purple', False: 'lightgrey'}, alpha=0.6)
    
    plt.title(title)
    plt.xlabel('Effect Size (Coefficient)')
    plt.ylabel('-log10(p-value)')
    plt.axhline(-np.log10(alpha), color='red', linestyle='--', alpha=0.5)
    
    safe_title = title.replace(" ", "_").replace("/", "-")
    fig_file = figures_dir / f'{project}.{modality}.{regression_type}_volcano.{safe_title}.png'
    plt.savefig(fig_file, dpi=300)
    plt.close()
    logger.info(f"Saved volcano plot to {fig_file}")


def main():
    args = parse_args()
    debug = args.debug

    # Configure logging
    log_filename = f"{args.modality}_{args.regression_type}_post_regression.log"
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
    results_dir = work_dir / "results"
    figures_dir = work_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    project = args.project
    modality = args.modality.lower() # ensure lowercase to match file patterns
    regression_type = args.regression_type
    
    # In the notebook, 'GEX' was used, but files seem to use lowercase 'rna'/'atac' based on previous scripts.
    # The user input args.modality is lowercase from choices. 
    # NOTE: The notebook used 'GEX'/'ATAC' uppercase for some vars but 'rna'/'atac' in filenames? 
    # Let's verify file pattern.
    # Notebook: results_file = f'{results_dir}/{project}.{modality}.{prefix_type}.{REGRESSION_TYPE}.age.csv'
    # Notebook Params: modality = 'GEX'
    # BUT prep_pb_data.py uses 'rna'/'atac'.
    # pseudobulk_regression.py uses 'rna'/'atac'.
    # So we should stick to lowercase 'rna'/'atac'.
    
    # 1. Aggregation
    logger.info("Aggregating regression results...")
    
    # Pattern to match regression result files
    # Format: {project}.{modality}.{cell_type}.{regression_type}.age.csv
    # e.g. aging_phase2.rna.Microglia.glm_tweedie.age.csv
    
    # We can glob them
    pattern = f"{project}.{modality}.*.{regression_type}.age.csv"
    result_files = list(results_dir.glob(pattern))
    
    if not result_files:
        logger.warning(f"No result files found matching {pattern} in {results_dir}")
        return

    dfs = []
    for file_path in result_files:
        try:
            df = read_csv(file_path)
            # cell_type is in the filename, but also typically in the dataframe column 'tissue'
            # if the regression script put it there.
            if 'tissue' not in df.columns:
                 # Extract from filename if missing: aging_phase2.rna.{cell_type}.glm_tweedie.age.csv
                 parts = file_path.name.split('.')
                 # Assuming project has 1 dot "aging_phase2" -> parts[0]
                 # Wait, project="aging_phase2" is one part.
                 # aging_phase2.rna.Microglia.glm_tweedie.age.csv
                 # parts: [0]=aging_phase2, [1]=rna, [2]=Microglia, [3]=glm_tweedie, [4]=age, [5]=csv
                 # This is risky if project has dots. 
                 # Better to trust the 'tissue' column from pseudobulk_regression.py which we know we added.
                 pass 
            
            dfs.append(df)
        except Exception as e:
            logger.error(f"Failed to read {file_path}: {e}")

    if not dfs:
        logger.error("No valid dataframes loaded.")
        return

    glm_results = concat(dfs, ignore_index=True)
    logger.info(f"Aggregated results shape: {glm_results.shape}")
    
    if debug:
        print(glm_results.head())

    # 2. FDR Correction
    logger.info("Computing BH FDR...")
    glm_results['p-value'] = glm_results['p-value'].fillna(1)
    glm_results = compute_bh_fdr(glm_results, verbose=True)

    # 3. RLM Effect Size Filtering
    if regression_type == 'rlm':
        logger.info(f"Applying RLM effect size filter (min_effect={args.min_rlm_effect})...")
        mask = glm_results.coef.abs() < args.min_rlm_effect
        n_filtered = mask.sum()
        glm_results.loc[mask, 'fdr_bh'] = 1.0
        logger.info(f"Set {n_filtered} results to FDR=1.0 due to small effect size.")
        
    # Count significant
    n_sig = (glm_results['fdr_bh'] < 0.05).sum()
    logger.info(f"Total significant features (FDR < 0.05): {n_sig}")

    # 4. Summary Counts
    total_counts = glm_results['tissue'].value_counts()
    sig_counts = glm_results.loc[glm_results['fdr_bh'] < 0.05]['tissue'].value_counts()
    
    summary_df = pd.DataFrame({'Total': total_counts, 'Significant': sig_counts})
    summary_df['Significant'] = summary_df['Significant'].fillna(0).astype(int)
    summary_df['Percentage'] = (summary_df['Significant'] / summary_df['Total'] * 100).round(2)
    summary_df = summary_df.sort_values('Percentage', ascending=False)
    
    logger.info("Significant counts and percentages per cell type:")
    print(summary_df)
    
    # 5. Save Results
    # Output filenames
    # notebook used: {project}.{modality}.{prefix_type}.{REGRESSION_TYPE}.age.csv
    # We don't have "prefix_type" (broad/specific) explicitly passed, usually implied by cell types?
    # Or maybe we should output one big file for everything. 
    # Notebook logic: 'curated_type' -> 'broad', 'cluster_name' -> 'specific'.
    # Our scripts seem to run per cell-type without knowing if it's broad or specific.
    # We will just name it "all_celltypes".
    
    results_file = results_dir / f'{project}.{modality}.all_celltypes.{regression_type}.age.csv'
    results_fdr_file = results_dir / f'{project}.{modality}.all_celltypes.{regression_type}_fdr.age.csv'
    
    logger.info(f"Saving full results to {results_file}")
    glm_results.to_csv(results_file, index=False)
    
    logger.info(f"Saving significant results to {results_fdr_file}")
    glm_results.loc[glm_results['fdr_bh'] < 0.05].to_csv(results_fdr_file, index=False)

    # 6. Volcano Plots
    logger.info("Generating volcano plots...")
    
    # Plot all results
    volcano_plot(glm_results, project, modality, regression_type, figures_dir, title='All Cell Types')
    
    # Plot per cell type
    cell_types = glm_results['tissue'].unique()
    for ct in cell_types:
        ct_results = glm_results.loc[glm_results.tissue == ct]
        if ct_results.shape[0] > 0:
            volcano_plot(ct_results, project, modality, regression_type, figures_dir, title=ct)


if __name__ == "__main__":
    main()
