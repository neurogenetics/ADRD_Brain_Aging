import argparse
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot cis correlation summary results by computing them from the full results."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--endo-modality", type=str, default="rna")
    parser.add_argument("--exog-modality", type=str, default="atac")
    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
        choices=["ols", "rlm", "glm", "glm_tweedie", "wls", "vwrlm"],
    )
    parser.add_argument("--fdr-threshold", type=float, default=0.05)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figs_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"

    log_filename = f"{logs_dir}/{args.project}_{args.endo_modality}_{args.exog_modality}_{args.regression_type}_cis_summary_plot.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )

    endo_results_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}.all_celltypes.{args.regression_type}_fdr_filtered.age.csv"
    )
    exog_results_file = (
        results_dir
        / f"{args.project}.all_celltypes.{args.exog_modality}.{args.regression_type}.age.csv"
    )
    cis_sig_file = (
        results_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.cis.fdr_filtered.csv"
    )

    for f in [endo_results_file, exog_results_file, cis_sig_file]:
        if not f.exists():
            logger.error(f"Required file not found: {f}")
            return

    logger.info(f"Loading endo results from {endo_results_file}")
    endo_results = pd.read_csv(endo_results_file)
    endo_results = endo_results[endo_results["fdr_bh"] < args.fdr_threshold]

    logger.info(f"Loading exog results from {exog_results_file}")
    exog_results = pd.read_csv(exog_results_file)
    # the exog file is not fdr filtered directly in the cis_correlation.py script usually
    # exog_results = exog_results[exog_results["fdr_bh"] < args.fdr_threshold]

    logger.info(f"Loading significant cis results from {cis_sig_file}")
    sig_results_df = pd.read_csv(cis_sig_file)

    cell_types = endo_results["tissue"].unique()
    
    plot_rows = []
    endo_upper = args.endo_modality.upper()
    exog_upper = args.exog_modality.upper()

    for cell_type in sorted(cell_types):
        ct_endo = endo_results[endo_results["tissue"] == cell_type]["feature"].nunique()
        ct_exog = exog_results[exog_results["tissue"] == cell_type]["feature"].nunique()

        ct_sig = sig_results_df[sig_results_df["tissue"] == cell_type]
        corr_endo = ct_sig["endo_feature"].nunique()
        corr_exog = ct_sig["exog_feature"].nunique()

        pct_endo = (corr_endo / ct_endo * 100) if ct_endo > 0 else 0
        pct_exog = (corr_exog / ct_exog * 100) if ct_exog > 0 else 0

        plot_rows.append(
            {
                "Cell Type": cell_type,
                "Percent Correlated": pct_endo,
                "Modality": endo_upper,
            }
        )
        plot_rows.append(
            {
                "Cell Type": cell_type,
                "Percent Correlated": pct_exog,
                "Modality": exog_upper,
            }
        )

    plot_df = pd.DataFrame(plot_rows)

    figs_dir.mkdir(parents=True, exist_ok=True)
    fig_file = (
        figs_dir
        / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.cis_summary.png"
    )

    logger.info("Generating plot...")
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))

    ax = sns.barplot(
        data=plot_df,
        x="Percent Correlated",
        y="Cell Type",
        hue="Modality",
        palette="colorblind",
    )
    ax.set_xlabel("Percent cis Correlated (%)")
    ax.set_ylabel("")
    plt.legend(title="Modality", loc="lower right")

    plt.tight_layout()
    plt.savefig(fig_file, dpi=300)
    plt.close()
    logger.info(f"Saved visualization to {fig_file}")


if __name__ == "__main__":
    main()