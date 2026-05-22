import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
from cnmf import cNMF

from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.linear_model import ElasticNetCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

from pseudobulk_convert import MODAL_TYPES_DICT

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Predict chronological age from age-associated latent factors."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument(
        "--density-threshold",
        type=float,
        default=0.1,
        help="Density threshold used in cNMF.",
    )
    parser.add_argument("--debug", action="store_true")
    parser.add_argument(
        "--cnmf-dir-name",
        type=str,
        default="cnmf",
        help="Name of the cnmf output directory within the latents path.",
    )
    parser.add_argument(
        "--p-threshold",
        type=float,
        default=0.05,
        help="P-value threshold for inclusion of factors.",
    )
    parser.add_argument(
        "--cv-folds",
        type=int,
        default=5,
        help="Number of folds for cross-validation.",
    )
    return parser.parse_args()


def load_sig_factors(results_dir, project, modality, p_threshold):
    fdr_file = (
        results_dir
        / "latents"
        / "cnmf"
        / "wls_pb_results"
        / f"{project}_combined_{modality}_pb_wls_fdr.csv"
    )

    if not fdr_file.exists():
        logger.warning(f"File not found: {fdr_file}")
        return pd.DataFrame()

    df = pd.read_csv(fdr_file)
    if "pval_age" not in df.columns:
        logger.warning(f"'pval_age' column missing in {fdr_file}")
        return pd.DataFrame()

    sig_df = df[df["pval_age"] <= p_threshold].copy()
    sig_df["modality"] = modality
    logger.info(
        f"Loaded {len(sig_df)} nominally significant factors (pval_age <= {p_threshold}) for {modality}."
    )
    return sig_df


def evaluate_model(model, X, y, cv, name, scope, modality):
    # Generate cross-validated predictions for scatter plots
    preds = cross_val_predict(model, X, y, cv=cv, n_jobs=-1)

    # For LOOCV, we calculate metrics globally across the out-of-fold predictions
    rmse = np.sqrt(mean_squared_error(y, preds))
    mae = mean_absolute_error(y, preds)
    r2 = r2_score(y, preds)

    return {
        "model": name,
        "scope": scope,
        "modality": modality,
        "rmse": rmse,
        "rmse_std": 0.0,
        "mae": mae,
        "r2": r2,
    }, preds


def main():
    args = parse_args()
    debug = args.debug

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    results_dir = work_dir / "results"
    logs_dir = work_dir / "logs"
    figs_dir = work_dir / "figures"

    out_dir = results_dir / "latents" / "cnmf" / "predictions"
    out_dir.mkdir(parents=True, exist_ok=True)

    logs_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    cnmf_dir = results_dir / "latents" / args.cnmf_dir_name

    log_filename = logs_dir / f"{args.project}_predict_age_from_latents.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
    )

    logger.info(f"Command line: {' '.join(sys.argv)}")

    # 1. Identify nominally significant factors
    rna_sig = load_sig_factors(results_dir, args.project, "rna", args.p_threshold)
    atac_sig = load_sig_factors(results_dir, args.project, "atac", args.p_threshold)

    sig_factors = pd.concat([rna_sig, atac_sig], ignore_index=True)
    if sig_factors.empty:
        logger.error("No significant factors found.")
        sys.exit(1)

    logger.info(f"Total significant factors to evaluate: {len(sig_factors)}")

    annot_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"
    logger.info(f"Loading annotated data obs from {annot_file}...")
    try:
        adata = sc.read_h5ad(annot_file, backed="r")
        # Load sample_id, age, cell_label, modality
        obs_df = adata.obs[["modality", "cell_label", "sample_id", "age"]].copy()
        obs_df["safe_cell_label"] = (
            obs_df["cell_label"].str.replace(" ", "_").str.replace("/", "-")
        )
        obs_df["age"] = pd.to_numeric(obs_df["age"], errors="coerce")
    except Exception as e:
        logger.error(f"Failed to load annotated anndata object: {e}")
        sys.exit(1)

    # 2. Extract Usages and build sample-level feature matrix
    factor_usages = {}
    grouped_factors = sig_factors.groupby(["cell_type", "modality", "k"])

    for (ct, mod, k), group in grouped_factors:
        run_name = f"{args.project}_{ct}_{mod}"
        logger.info(f"Loading cNMF usages for {run_name} K={k}")

        cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=run_name)
        try:
            usage, *_ = cnmf_obj.load_results(
                K=k, density_threshold=args.density_threshold
            )
        except Exception as e:
            logger.error(
                f"Failed to load cNMF results for {run_name} K={k}. Error: {e}"
            )
            continue

        modality_list = MODAL_TYPES_DICT.get(mod, [mod])
        ct_obs = obs_df[
            (obs_df.modality.isin(modality_list)) & (obs_df["safe_cell_label"] == ct)
        ].copy()
        ct_obs["sample_id"] = ct_obs["sample_id"].astype(str)

        usage = usage.reindex(ct_obs.index).dropna()
        if len(usage) == 0:
            continue

        usage_with_sample = usage.merge(
            ct_obs[["sample_id"]], left_index=True, right_index=True, how="left"
        )
        usage_with_sample = usage_with_sample.dropna(subset=["sample_id"])

        # Pseudobulk to sample level
        agg_usage = usage_with_sample.groupby("sample_id").mean()

        for factor_idx in group["factor"].values:
            factor_str = str(factor_idx)
            factor_int = int(factor_idx)

            if factor_str in agg_usage.columns:
                col_name = f"{ct}|{mod.upper()}|K{k}|F{factor_str}"
                factor_usages[col_name] = agg_usage[factor_str]
            elif factor_int in agg_usage.columns:
                col_name = f"{ct}|{mod.upper()}|K{k}|F{factor_idx}"
                factor_usages[col_name] = agg_usage[factor_int]

    if not factor_usages:
        logger.error("No usage data could be extracted.")
        sys.exit(1)

    X_full = pd.DataFrame(factor_usages).fillna(0)

    # Get age target vector
    sample_ages = (
        obs_df.dropna(subset=["sample_id", "age"]).groupby("sample_id")["age"].first()
    )

    # Align X and y
    common_samples = X_full.index.intersection(sample_ages.index)
    X_full = X_full.loc[common_samples]
    y_full = sample_ages.loc[common_samples]

    logger.info(f"Final feature matrix shape: {X_full.shape}")

    # Set up models
    models = {
        "ElasticNet": Pipeline(
            [
                ("scaler", StandardScaler()),
                (
                    "enet",
                    ElasticNetCV(
                        cv=3,
                        l1_ratio=[0.1, 0.5, 0.7, 0.9, 0.95, 0.99, 1],
                        random_state=42,
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
        "RandomForest": RandomForestRegressor(
            n_estimators=100, random_state=42, n_jobs=-1
        ),
    }

    cv = LeaveOneOut()

    results = []
    predictions = {}  # Store cross-val predictions for the best models

    # --- A. INDIVIDUAL FACTOR EVALUATION ---
    logger.info("Evaluating individual factors...")
    for col in X_full.columns:
        ct, mod, _, _ = col.split("|")
        X_ind = X_full[[col]]
        for model_name, model in models.items():
            res, preds = evaluate_model(
                model, X_ind, y_full, cv, model_name, f"Individual ({col})", mod
            )
            results.append(res)

    # --- B. CELL-TYPE GROUP EVALUATION ---
    logger.info("Evaluating cell-type groups...")
    cell_types = list(set([col.split("|")[0] for col in X_full.columns]))
    for ct in cell_types:
        for mod in ["RNA", "ATAC"]:
            cols = [c for c in X_full.columns if c.startswith(f"{ct}|{mod}")]
            if len(cols) > 0:
                X_ct = X_full[cols]
                for model_name, model in models.items():
                    res, preds = evaluate_model(
                        model,
                        X_ct,
                        y_full,
                        cv,
                        model_name,
                        f"Cell-Type ({ct}) - {mod}",
                        mod,
                    )
                    results.append(res)
                    predictions[f"Cell-Type ({ct}) - {mod} - {model_name}"] = preds

    # --- C. GLOBAL EVALUATION ---
    logger.info("Evaluating global models...")
    for mod in ["RNA", "ATAC"]:
        cols = [c for c in X_full.columns if f"|{mod}|" in c]
        if len(cols) > 0:
            X_mod = X_full[cols]
            for model_name, model in models.items():
                res, preds = evaluate_model(
                    model,
                    X_mod,
                    y_full,
                    cv,
                    model_name,
                    f"Global (All Cell Types) - {mod}",
                    mod,
                )
                results.append(res)
                predictions[f"Global (All Cell Types) - {mod} - {model_name}"] = preds

                # Fit on full data to get feature importances
                model.fit(X_mod, y_full)
                if model_name == "ElasticNet":
                    importances = pd.DataFrame(
                        {
                            "Factor": cols,
                            "Importance (Coef)": model.named_steps["enet"].coef_,
                        }
                    )
                else:
                    importances = pd.DataFrame(
                        {
                            "Factor": cols,
                            "Importance (Feature Imp)": model.feature_importances_,
                        }
                    )

                imp_out = (
                    out_dir
                    / f"{args.project}_age_prediction_importances_{mod}_{model_name}.csv"
                )
                importances.sort_values(
                    by=importances.columns[1], key=abs, ascending=False
                ).to_csv(imp_out, index=False)

    # --- COMBINED CROSS-MODALITY GLOBAL EVALUATION ---
    if (
        sum(1 for c in X_full.columns if "|RNA|" in c) > 0
        and sum(1 for c in X_full.columns if "|ATAC|" in c) > 0
    ):
        logger.info("Evaluating combined cross-modality global models...")
        for model_name, model in models.items():
            res, preds = evaluate_model(
                model,
                X_full,
                y_full,
                cv,
                model_name,
                "Global Cross-Modality",
                "RNA+ATAC",
            )
            results.append(res)
            predictions[f"Global Cross-Modality - {model_name}"] = preds

            model.fit(X_full, y_full)
            if model_name == "ElasticNet":
                importances = pd.DataFrame(
                    {
                        "Factor": X_full.columns,
                        "Importance (Coef)": model.named_steps["enet"].coef_,
                    }
                )
            else:
                importances = pd.DataFrame(
                    {
                        "Factor": X_full.columns,
                        "Importance (Feature Imp)": model.feature_importances_,
                    }
                )

            imp_out = (
                out_dir
                / f"{args.project}_age_prediction_importances_Multiome_{model_name}.csv"
            )
            importances.sort_values(
                by=importances.columns[1], key=abs, ascending=False
            ).to_csv(imp_out, index=False)

    # --- SAVE RESULTS & PLOT ---
    results_df = pd.DataFrame(results)
    results_out = out_dir / f"{args.project}_age_prediction_results.csv"
    results_df.sort_values("rmse").to_csv(results_out, index=False)
    logger.info(f"Saved prediction results to {results_out}")

    # Plot top 15 models by RMSE
    try:
        # Group by scope to get the top 15 scopes with the best minimum RMSE
        top_scopes = (
            results_df.groupby("scope")["rmse"].min().sort_values().head(15).index
        )
        top_models = results_df[results_df["scope"].isin(top_scopes)]

        plt.figure(figsize=(12, 8))
        sns.barplot(
            data=top_models,
            x="rmse",
            y="scope",
            hue="model",
            dodge=True,
            palette="colorblind",
            order=top_scopes,
        )
        plt.title("Top 15 Age Prediction Scopes (by Cross-Validated RMSE)")
        plt.xlabel("RMSE (Years)")
        plt.ylabel("Model Scope")
        plt.tight_layout()
        plt.savefig(figs_dir / f"{args.project}_top_age_prediction_models.png", dpi=300)
        plt.close()
    except Exception as e:
        logger.error(f"Failed to plot summary: {e}")

    # Plot best actual vs predicted scatter
    try:
        best_row = results_df.sort_values("rmse").iloc[0]
        best_name = (
            f"{best_row['scope']} - {best_row['modality']} - {best_row['model']}"
        )
        # Fallback naming for cross-modality
        if best_row["scope"] == "Global Cross-Modality":
            best_name = f"Global Cross-Modality - {best_row['model']}"

        if best_name in predictions:
            best_preds = predictions[best_name]
            plt.figure(figsize=(8, 8))
            sns.scatterplot(x=y_full, y=best_preds, alpha=0.7)

            # Line of perfect prediction
            min_val = min(y_full.min(), min(best_preds))
            max_val = max(y_full.max(), max(best_preds))
            plt.plot([min_val, max_val], [min_val, max_val], "r--")

            plt.title(
                f"Best Age Prediction: {best_row['scope']}\nModel: {best_row['model']} | RMSE: {best_row['rmse']:.2f}"
            )
            plt.xlabel("Actual Age")
            plt.ylabel("Predicted Age")
            plt.tight_layout()
            plt.savefig(
                figs_dir / f"{args.project}_best_age_prediction_scatter.png", dpi=300
            )
            plt.close()
    except Exception as e:
        logger.error(f"Failed to plot scatter: {e}")

    logger.info("Age prediction analysis complete.")


if __name__ == "__main__":
    main()
