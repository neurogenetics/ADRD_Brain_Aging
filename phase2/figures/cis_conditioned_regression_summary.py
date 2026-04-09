import argparse
import logging
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Summarize and visualize the cis conditioned regression analysis."
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
        "--alpha", 
        type=float, 
        default=0.05, 
        help="Alpha threshold for mediation (uncorrected p-value of age exposure)"
    )
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figures_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"

    figures_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_filename = logs_dir / f"{args.project}_{args.endo_modality}_{args.exog_modality}_{args.regression_type}_conditioned_summary.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    # Input files
    endo_fdr_file = results_dir / f"{args.project}.{args.endo_modality}.all_celltypes.{args.regression_type}_fdr_filtered.age.csv"
    cis_fdr_file = results_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.cis.fdr_filtered.csv"
    cond_file = results_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.conditioned.csv"

    for f in [endo_fdr_file, cis_fdr_file, cond_file]:
        if not f.exists():
            logger.error(f"Required file not found: {f}")
            return

    logger.info("Loading required result files...")
    endo_df = pd.read_csv(endo_fdr_file)
    cis_df = pd.read_csv(cis_fdr_file)
    cond_df = pd.read_csv(cond_file)

    logger.info(f"Loaded {len(endo_df)} endo age-associated features.")
    logger.info(f"Loaded {len(cis_df)} significant cis-correlated pairs.")
    logger.info(f"Loaded {len(cond_df)} conditioned pairs.")

    summary_rows = []
    
    logger.info("Computing summary metrics per gene and cell type...")
    for tissue, group in endo_df.groupby("tissue"):
        tissue_genes = group["feature"].unique()
        
        tissue_cis = cis_df[cis_df["tissue"] == tissue]
        tissue_cond = cond_df[cond_df["tissue"] == tissue]
        
        for gene in tissue_genes:
            cis_cor_peaks = tissue_cis[tissue_cis["endo_feature"] == gene].shape[0]
            
            this_cond = tissue_cond[tissue_cond["endo_feature"] == gene]
            cis_cor_age_peaks = this_cond.shape[0]
            
            # Mediated if conditioned exposure p-value > alpha (loss of significance for age)
            mediated = this_cond[this_cond["exposure_pval"] > args.alpha]
            mediating_peak_count = mediated.shape[0]
            
            summary_rows.append({
                "feature": gene,
                "tissue": tissue,
                "cis_cor_peaks": cis_cor_peaks,
                "cis_cor_age_peaks": cis_cor_age_peaks,
                "mediating_peak_count": mediating_peak_count
            })
            
    summary_df = pd.DataFrame(summary_rows)

    # Compute proportions table
    this_list = []
    for tissue in summary_df["tissue"].unique():
        tissue_summary = summary_df[summary_df["tissue"] == tissue]
        cell_type_cnt = tissue_summary["feature"].nunique()
        mediated_cnt = tissue_summary[tissue_summary["mediating_peak_count"] > 0].shape[0]
        mediated_percent = (mediated_cnt / cell_type_cnt * 100) if cell_type_cnt > 0 else 0
        this_list.append({
            "tissue": tissue,
            "count": cell_type_cnt,
            "percent": mediated_percent,
            "pairwise_cnt": mediated_cnt
        })
    mediated_proportions = pd.DataFrame(this_list)

    # Save summary table
    out_summary_file = figures_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.conditioned.age.summary.csv"
    mediated_proportions.to_csv(out_summary_file, index=False)
    logger.info(f"Saved summary proportions table to {out_summary_file}")

    # File prefixes for figures
    fig_scatter = figures_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.conditioned.summary"
    fig_bar = figures_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.conditioned.summary_bar"
    fig_dist = figures_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.conditioned.summary_dist"
    fig_cnt = figures_dir / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.conditioned.summary_cnt"

    sns.set_theme(style='whitegrid')

    # 1. Scatter Plot
    logger.info("Generating scatter plot...")
    plt.figure(figsize=(11, 11), dpi=100)
    sns.scatterplot(
        data=summary_df.sample(frac=1, random_state=42), 
        x='cis_cor_peaks', 
        y='mediating_peak_count', 
        hue='tissue', 
        size='cis_cor_age_peaks', 
        palette='bright'
    )
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
    plt.tight_layout()
    plt.savefig(f"{fig_scatter}.png", bbox_inches='tight')
    plt.savefig(f"{fig_scatter}.svg", bbox_inches='tight')
    plt.close()

    # 2. Bar Plot
    logger.info("Generating bar plot...")
    plt.figure(figsize=(15, 11), dpi=100)
    sns.barplot(
        data=mediated_proportions.sort_values('percent', ascending=False),
        x='tissue', 
        y='percent', 
        color='purple'
    )
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.title('% of age associated genes that are mediated by a cis correlated age associated ATAC peak', fontsize='large')
    plt.xlabel('Cell types')
    plt.ylabel('%')
    plt.savefig(f"{fig_bar}.png", bbox_inches='tight')
    plt.savefig(f"{fig_bar}.svg", bbox_inches='tight')
    plt.close()

    # Distance calculations
    distances = []
    sig_cond = cond_df[cond_df["exposure_pval"] > args.alpha]
    
    for row in summary_df[summary_df["mediating_peak_count"] > 1].itertuples():
        this_mediated = sig_cond[(sig_cond["endo_feature"] == row.feature) & (sig_cond["tissue"] == row.tissue)]
        if len(this_mediated) > 1:
            exog_features = this_mediated["exog_feature"]
            starts = []
            ends = []
            for peak in exog_features:
                try:
                    parts = peak.split(":")[1].split("-")
                    starts.append(int(parts[0]))
                    ends.append(int(parts[1]))
                except Exception:
                    continue
            
            if len(starts) > 1:
                starts = np.array(starts)
                ends = np.array(ends)
                midpoints = (starts + ends) / 2
                midpoints.sort()
                mean_dist = int(np.diff(midpoints).mean() / 1000)
                distances.append({
                    "feature": row.feature,
                    "tissue": row.tissue,
                    "mean_dist": mean_dist
                })
            
    distances_df = pd.DataFrame(distances)

    label_order = mediated_proportions.sort_values('percent', ascending=False)["tissue"].values

    # 3. Distance Boxen Plot
    if not distances_df.empty:
        logger.info("Generating distance boxen plot...")
        distances_df = distances_df.sort_values('mean_dist', ascending=False)
        
        plt.figure(figsize=(9, 9), dpi=100)
        try:
            sns.boxenplot(
                x='tissue', y='mean_dist', data=distances_df, 
                color='purple', order=label_order, width_method='exponential', k_depth='trustworthy'
            )
        except TypeError:
            sns.boxenplot(
                x='tissue', y='mean_dist', data=distances_df, 
                color='purple', order=label_order
            )
            
        sns.stripplot(
            x='tissue', y='mean_dist', data=distances_df, 
            alpha=0.75, jitter=True, color='darkgrey', order=label_order
        )
        plt.xticks(rotation=90)
        plt.title('Mean distances between cis proximal peaks that mediate age effects', fontsize='large')
        plt.xlabel('Cell types')
        plt.ylabel('Mean distance (Kb)')
        plt.tight_layout()
        plt.savefig(f"{fig_dist}.png", bbox_inches='tight')
        plt.savefig(f"{fig_dist}.svg", bbox_inches='tight')
        plt.close()
    else:
        logger.info("No distance data to plot (no genes with >1 mediating peak).")

    # 4. Peak Count Boxen Plot
    mediating_df = summary_df[summary_df["mediating_peak_count"] > 0].sort_values('mediating_peak_count', ascending=False)
    if not mediating_df.empty:
        logger.info("Generating peak count boxen plot...")
        plt.figure(figsize=(9, 9), dpi=100)
        try:
            sns.boxenplot(
                x='tissue', y='mediating_peak_count', data=mediating_df, 
                color='purple', order=label_order, width_method='exponential', k_depth='trustworthy'
            )
        except TypeError:
            sns.boxenplot(
                x='tissue', y='mediating_peak_count', data=mediating_df, 
                color='purple', order=label_order
            )
            
        sns.stripplot(
            x='tissue', y='mediating_peak_count', data=mediating_df, 
            alpha=0.75, jitter=True, color='darkgrey', order=label_order
        )
        plt.ylim(bottom=0)
        plt.xticks(rotation=90)
        plt.title('Number of cis proximal peaks that mediate age effects', fontsize='large')
        plt.xlabel('Cell types')
        plt.ylabel('Number mediating cis proximal ATAC')
        plt.tight_layout()
        plt.savefig(f"{fig_cnt}.png", bbox_inches='tight')
        plt.savefig(f"{fig_cnt}.svg", bbox_inches='tight')
        plt.close()

    logger.info("Summary and plotting completed.")


if __name__ == "__main__":
    main()
