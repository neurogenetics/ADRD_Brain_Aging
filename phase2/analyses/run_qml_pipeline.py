import argparse
import logging
from pathlib import Path
import sys

import pandas as pd
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set PYTORCH_CUDA_ALLOC_CONF programmatically before anything is allocated
import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

# Add src to path so we can import qml_pipeline
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from qml_pipeline import (
    ClassicalWalk, 
    QuantumWalk, 
    CompletePipeline, 
    train_pipeline,
    extract_key_features_kneed,
    compare_walk_features
)

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"

def parse_args():
    parser = argparse.ArgumentParser(description="Run the Quantum Bioinformatics Pipeline (CTQW/CTRW) on cNMF results.")
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--cell-type", type=str, required=True, help="Target cell type.")
    parser.add_argument("--modality", type=str, default="rna", choices=["rna", "atac"])
    parser.add_argument("--cnmf-dir-name", type=str, default="cnmf", help="Name of the cnmf output directory within the latents path.")
    parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs.")
    parser.add_argument("--lr", type=float, default=0.01, help="Learning rate.")
    parser.add_argument("--seed", type=int, default=42, help="Random seed for reproducibility.")
    parser.add_argument("--max-nodes", type=int, default=16384, help="Maximum number of features to process to prevent OOM errors (VRAM protection).")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()

def get_target_age_for_cell_type(info_dir: Path, project: str, cell_type: str, modality: str) -> float:
    """Gets the mean age of the donors for this cell type to use as a dummy/baseline target for the pipeline run."""
    safe_ct = cell_type.replace(" ", "_").replace("/", "-")
    covars_file = info_dir / f"{project}.{safe_ct}.{modality}.final_covariates.csv"
    if not covars_file.exists():
        logger.warning(f"Covariates file not found: {covars_file}. Using default age 65.0")
        return 65.0
    
    df = pd.read_csv(covars_file, index_col=0)
    if "age" in df.columns:
        return df["age"].mean()
    return 65.0

def save_full_walk_results(probs: torch.Tensor, feature_names: list, key_features: list, base_out_file: Path):
    """
    Saves the full list of features sorted by their final walk probability,
    and additionally saves a subset file containing only the extracted 'key_features'.
    """
    probs_np = probs.cpu().numpy()
    
    # Create DataFrame
    df = pd.DataFrame({
        "feature": feature_names,
        "probability": probs_np
    })
    
    # Sort descending by probability
    df = df.sort_values(by="probability", ascending=False).reset_index(drop=True)
    
    # Save the full distribution
    full_out_file = base_out_file.parent / f"{base_out_file.stem}_full_probs.csv"
    df.to_csv(full_out_file, index=False)
    logger.info(f"Saved full sorted walk probabilities to {full_out_file}")
    
    # Subset to only the 'key' features extracted by kneed
    key_df = df[df["feature"].isin(key_features)].reset_index(drop=True)
    
    # Save the key feature subset
    key_out_file = base_out_file.parent / f"{base_out_file.stem}_key_features.csv"
    key_df.to_csv(key_out_file, index=False)
    logger.info(f"Saved key feature subset ({len(key_df)} features) to {key_out_file}")

def generate_probability_scatter_plot(c_probs: torch.Tensor, q_probs: torch.Tensor, feature_names: list, classical_keys: list, quantum_keys: list, out_file: Path, title: str):
    """Generates a scatter plot comparing Classical vs Quantum probabilities for all features."""
    c_probs_np = c_probs.cpu().numpy()
    q_probs_np = q_probs.cpu().numpy()
    
    df = pd.DataFrame({
        "feature": feature_names,
        "Classical Probability": c_probs_np,
        "Quantum Probability": q_probs_np
    })
    
    c_set = set(classical_keys)
    q_set = set(quantum_keys)
    
    def get_category(feature):
        if feature in c_set and feature in q_set:
            return "Shared"
        elif feature in c_set:
            return "Classical Specific"
        elif feature in q_set:
            return "Quantum Specific"
        return "Background"
    
    df["Category"] = df["feature"].apply(get_category)
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(
        data=df, 
        x="Classical Probability", 
        y="Quantum Probability", 
        hue="Category",
        palette={"Classical Specific": "red", "Quantum Specific": "blue", "Shared": "purple", "Background": "lightgrey"},
        alpha=0.6,
        s=20
    )
    
    # Draw y=x diagonal
    max_val = max(df["Classical Probability"].max(), df["Quantum Probability"].max())
    plt.plot([0, max_val], [0, max_val], color='black', linestyle='--', alpha=0.5)
    
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Saved Probability Scatter Plot to {out_file}")

def main():
    args = parse_args()
    
    import random
    seed = args.seed
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    info_dir = work_dir / "sample_info"
    figs_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"
    
    figs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    
    log_filename = logs_dir / f"{args.project}_{args.modality}_{safe_ct}_qml_pipeline.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    
    logger.info(f"Command line: {' '.join(sys.argv)}")
    
    cnmf_dir = results_dir / "latents" / args.cnmf_dir_name
    out_dir = cnmf_dir / "qml_results"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Determine selected K from post_cnmf_latent_regressions.py outputs
    regression_type = "pb_wls" # Defaulting to pseudobulk WLS
    results_target_dir = cnmf_dir / "wls_pb_results"
    reg_file = results_target_dir / f"{args.project}_{safe_ct}_{args.modality}_{regression_type}.csv"
    
    if not reg_file.exists():
        logger.error(f"Regression results file not found: {reg_file}. Run cnmf_latent_regressions.py first.")
        sys.exit(1)
        
    reg_df = pd.read_csv(reg_file)
    if reg_df.empty:
        logger.error(f"Regression results file is empty: {reg_file}")
        sys.exit(1)
        
    selected_k = int(reg_df["k"].iloc[0])
    logger.info(f"Loaded selected K={selected_k} from regression results.")
    
    # 2. Load cNMF Spectra Matrix (H matrix)
    run_name = f"{args.project}_{safe_ct}_{args.modality}"
    run_dir = cnmf_dir / run_name
    
    score_files = list(run_dir.glob(f"*spectra_score.k_{selected_k}.dt_*.txt"))
    if not score_files:
        logger.error(f"No spectra score file found for K={selected_k} in {run_dir}.")
        sys.exit(1)
        
    score_file = score_files[0]
    logger.info(f"Loading cNMF spectra scores from: {score_file}")
    h_df = pd.read_csv(score_file, sep='\t', index_col=0)
    
    # 2.5 VRAM Protection: Dynamic Dimensionality Reduction
    if h_df.shape[1] > args.max_nodes:
        logger.warning(f"Feature count ({h_df.shape[1]}) exceeds safety limit ({args.max_nodes}). Downsampling to protect VRAM...")
        # Calculate variance for each feature across the K latent factors
        feature_variances = h_df.var(axis=0)
        # Select the top N features with the highest variance
        top_features = feature_variances.nlargest(args.max_nodes).index
        # Subset the H matrix
        h_df = h_df[top_features]
        logger.info(f"Dynamically downsampled to the {args.max_nodes} most variable features.")
    
    # H matrix shape is (factors, features)
    feature_names = h_df.columns.tolist()
    logger.info(f"Loaded H matrix with shape: {h_df.shape} (factors: {h_df.shape[0]}, features: {h_df.shape[1]})")
    
    # Convert to PyTorch Tensor
    h_tensor = torch.tensor(h_df.values, dtype=torch.float32)
    h_matrices = {f"{safe_ct}_{args.modality}": h_tensor}
    
    # 3. Get Target Age
    target_age = get_target_age_for_cell_type(info_dir, args.project, args.cell_type, args.modality)
    logger.info(f"Using target mean age: {target_age:.2f} for training.")
    
    device = "cuda" if torch.cuda.is_available() else "cpu"
    logger.info(f"Using device: {device}")
    
    # 4. Train Classical Pipeline
    logger.info("--- Initializing Classical CTRW Pipeline ---")
    classical_walk = ClassicalWalk(h_matrices, device=device)
    classical_pipeline = CompletePipeline(classical_walk)
    
    logger.info("Training Classical Pipeline...")
    classical_results = train_pipeline(
        classical_pipeline, 
        target_age=target_age, 
        num_epochs=args.epochs, 
        lr=args.lr
    )
    
    with torch.no_grad():
        classical_probs = classical_pipeline.walk_module()
        
    # Extract Key features and Save full probabilities
    classical_features = extract_key_features_kneed(classical_probs, feature_names)
    logger.info(f"Extracted {len(classical_features)} key features from Classical Walk.")
    
    c_out_file = out_dir / f"{args.project}_{safe_ct}_{args.modality}_classical"
    save_full_walk_results(classical_probs, feature_names, classical_features, c_out_file)
    
    # Explicitly clear classical models before instantiating quantum to save VRAM
    del classical_walk
    del classical_pipeline
    torch.cuda.empty_cache()
    
    # 5. Train Quantum Pipeline
    logger.info("--- Initializing Quantum CTQW Pipeline ---")
    quantum_walk = QuantumWalk(h_matrices, device=device)
    quantum_pipeline = CompletePipeline(quantum_walk)
    
    logger.info("Training Quantum Pipeline...")
    quantum_results = train_pipeline(
        quantum_pipeline, 
        target_age=target_age, 
        num_epochs=args.epochs, 
        lr=args.lr
    )
    
    with torch.no_grad():
        quantum_probs = quantum_pipeline.walk_module()
        
    # Extract Key features and Save full probabilities
    quantum_features = extract_key_features_kneed(quantum_probs, feature_names)
    logger.info(f"Extracted {len(quantum_features)} key features from Quantum Walk.")
    
    q_out_file = out_dir / f"{args.project}_{safe_ct}_{args.modality}_quantum"
    save_full_walk_results(quantum_probs, feature_names, quantum_features, q_out_file)
    
    # 6. Compare and Save
    logger.info("--- Comparing Walk Features ---")
    comparison = compare_walk_features(classical_features, quantum_features)
    
    logger.info(f"Jaccard Similarity: {comparison['jaccard_similarity']:.4f}")
    logger.info(f"Classical Specific Features: {len(comparison['classical_specific'])}")
    logger.info(f"Quantum Specific Features: {len(comparison['quantum_specific'])}")
    logger.info(f"Shared Features: {comparison['num_shared']}")
    
    comp_df = pd.DataFrame([comparison])
    comp_df["classical_specific"] = ",".join(comparison["classical_specific"])
    comp_df["quantum_specific"] = ",".join(comparison["quantum_specific"])
    
    out_file = out_dir / f"{args.project}_{safe_ct}_{args.modality}_qml_comparison.csv"
    comp_df.to_csv(out_file, index=False)
    logger.info(f"Saved QML comparison results to {out_file}")
    
    # 7. Generate Scatter Plot
    scatter_out = figs_dir / f"{args.project}_{safe_ct}_{args.modality}_qml_scatter.png"
    title = f"Quantum vs Classical Walk Probabilities\n{args.cell_type} ({args.modality.upper()})"
    generate_probability_scatter_plot(classical_probs, quantum_probs, feature_names, classical_features, quantum_features, scatter_out, title)

if __name__ == "__main__":
    main()
