import argparse
import logging
import sys
import os
from math import sqrt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from pandas import read_csv, read_parquet, DataFrame
from sklearn.linear_model import LogisticRegression
from sklearn.decomposition import PCA, FastICA, NMF
from matplotlib.pyplot import rc_context
from sklearn.metrics import r2_score, mean_squared_error, accuracy_score


def plot_pair(
    this_df: DataFrame,
    first: str,
    second: str,
    out_file: str,
    hue_cov=None,
    style_cov=None,
    size_cov=None,
):
    with rc_context({"figure.figsize": (8, 8), "figure.dpi": 100}):
        sns.set_style("whitegrid")
        sns_plot = sns.scatterplot(
            x=first,
            y=second,
            hue=hue_cov,
            style=style_cov,
            size=size_cov,
            data=this_df,
            palette="viridis",
        )
        plt.xlabel(first)
        plt.ylabel(second)
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, prop={"size": 10})
        plt.tight_layout()
        plt.savefig(out_file)
        plt.close()


def generate_selected_model(
    n_comps: int, data_df: DataFrame, model_type: str = "PCA"
) -> tuple:
    if model_type == "PCA":
        model = PCA(n_components=n_comps, random_state=42)
    elif model_type == "NMF":
        model = NMF(n_components=n_comps, init="random", random_state=42, max_iter=500)
    elif model_type == "ICA":
        model = FastICA(n_components=n_comps, random_state=42)
    
    components = model.fit_transform(data_df)
    recon_input = model.inverse_transform(components)
    r2 = r2_score(y_true=data_df, y_pred=recon_input)
    rmse = sqrt(mean_squared_error(data_df, recon_input))
    logging.debug(
        f"{model_type} with {n_comps} components accuracy is {r2:.4f}, RMSE is {rmse:.4f}"
    )
    ret_df = DataFrame(data=components, index=data_df.index).round(4)
    ret_df.columns = [f"{model_type}_{i}" for i in range(n_comps)]
    return model, ret_df, r2, rmse


def find_sex_mismatches(
    lat_df: DataFrame, target_df: DataFrame, target_col_name: str = "msex"
) -> tuple:
    comb_df = lat_df.join(target_df[[target_col_name]], how="left")
    # clean out nans
    comb_df = comb_df.dropna(subset=[target_col_name])
    if comb_df.empty:
        logging.warning(
            "No valid rows for sex mismatch prediction after dropping NaNs."
        )
        return [], comb_df

    logistic_model = LogisticRegression()
    logistic_model.fit(comb_df[lat_df.columns], comb_df[target_col_name])

    comb_df["pred"] = logistic_model.predict(comb_df[lat_df.columns])
    
    # Calculate probabilities
    probs = logistic_model.predict_proba(comb_df[lat_df.columns])
    for i, cls_name in enumerate(logistic_model.classes_):
        comb_df[f"prob_{cls_name}"] = probs[:, i]
    
    acc = accuracy_score(comb_df[target_col_name], comb_df["pred"])
    logging.info(
        f"logistic {target_col_name} predictor has {acc:.2%} accuracy"
    )
    
    sex_mismatch = comb_df.loc[
        comb_df[target_col_name] != comb_df["pred"]
    ]
    ret_list = list(sex_mismatch.index.values)
    return ret_list, comb_df


def main():
    parser = argparse.ArgumentParser(description="Run RNA sex check on sample data.")
    parser.add_argument("--rna_matrix", required=True, help="Input RNA matrix parquet file")
    parser.add_argument("--donor_info", required=True, help="Input donor information CSV file")
    parser.add_argument("--out_dir", required=True, help="Output directory for generated files")
    parser.add_argument("--sample_id_col", default=None, help="Column name in donor info CSV for sample IDs (default: use CSV index)")
    parser.add_argument("--sex_col", default="sex", help="Column name in donor info CSV for sex")
    parser.add_argument("--project", default="", help="Project prefix for output files")
    
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    os.makedirs(args.out_dir, exist_ok=True)

    logging.info(f"Loading RNA matrix from {args.rna_matrix}...")
    quants_df = read_parquet(args.rna_matrix)

    logging.info(f"Loading donor info from {args.donor_info}...")
    if args.sample_id_col:
        info_df = read_csv(args.donor_info)
        if args.sample_id_col not in info_df.columns:
            logging.error(f"Sample ID column '{args.sample_id_col}' not found in donor info CSV.")
            sys.exit(1)
        info_df = info_df.set_index(args.sample_id_col)
    else:
        info_df = read_csv(args.donor_info, index_col=0)

    sex_specific_genes = ["XIST", "RPS4Y1", "RPS4Y2", "KDM5D", "UTY", "DDX3Y", "USP9Y"]
    sex_genes_present = list(set(sex_specific_genes) & set(quants_df.columns))
    
    if not sex_genes_present:
        logging.warning("No sex-specific genes found in the RNA matrix.")
        sys.exit(0)

    quants_sex_df = quants_df[sex_genes_present]
    
    # Fill any NaNs in the RNA matrix with 0 before running models
    if quants_sex_df.isnull().values.any():
        logging.warning("NaNs detected in RNA matrix. Filling NaNs with 0.")
        quants_sex_df = quants_sex_df.fillna(0)
    
    if quants_sex_df.shape[1] > 0 and args.sex_col in info_df.columns:
        logging.info("Generating PCA model on sex-specific genes...")
        _, sex_pca_df, _, _ = generate_selected_model(2, quants_sex_df, "PCA")
        
        prefix = f"{args.project}_" if args.project else ""
        
        plot_out = os.path.join(args.out_dir, f"{prefix}sex_pca.png")
        logging.info(f"Plotting PCA results to {plot_out}...")
        
        plot_df = sex_pca_df.join(info_df[[args.sex_col]], how="left")
        plot_pair(
            plot_df,
            "PCA_0",
            "PCA_1",
            plot_out,
            hue_cov=args.sex_col,
        )
        
        logging.info("Predicting sex mismatches...")
        ids_sex_mismatch, predictions_df = find_sex_mismatches(sex_pca_df, info_df, target_col_name=args.sex_col)
        logging.info(f"Samples with predicted sex mismatch: {ids_sex_mismatch}")
        
        mismatch_out = os.path.join(args.out_dir, f"{prefix}mismatched_samples.txt")
        with open(mismatch_out, "w") as f:
            for sample_id in ids_sex_mismatch:
                f.write(f"{sample_id}\n")
        logging.info(f"Wrote mismatched sample IDs to {mismatch_out}")
        
        if not predictions_df.empty:
            predictions_out = os.path.join(args.out_dir, f"{prefix}sex_predictions.csv")
            predictions_df.to_csv(predictions_out)
            logging.info(f"Wrote complete sex predictions and probabilities to {predictions_out}")

    else:
        if args.sex_col not in info_df.columns:
            logging.error(f"Sex column '{args.sex_col}' not found in donor info CSV.")
        else:
            logging.error("Failed to extract sex-specific features.")
        sys.exit(1)


if __name__ == "__main__":
    main()
