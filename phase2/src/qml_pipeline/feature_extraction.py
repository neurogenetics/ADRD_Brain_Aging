import torch
import numpy as np
from kneed import KneeLocator
from typing import List

def extract_key_features_kneed(
    probabilities: torch.Tensor, 
    feature_names: List[str],
    S: float = 1.0
) -> List[str]:
    """
    Extracts the most influential features (genes/peaks) based on the point of 
    maximum curvature in the sorted probability distribution.
    
    Args:
        probabilities: 1D Tensor of state probabilities from the walk.
        feature_names: List of strings corresponding to the nodes.
        S: Sensitivity parameter for kneed.
    
    Returns:
        List of feature names that are highly influential.
    """
    if len(probabilities.shape) > 1:
        probabilities = probabilities.squeeze()
        
    probs_np = probabilities.detach().cpu().numpy()
    
    # Sort probabilities descending
    sorted_indices = np.argsort(probs_np)[::-1]
    sorted_probs = probs_np[sorted_indices]
    
    # X-axis is just the rank
    x = np.arange(len(sorted_probs))
    
    # Find the knee (elbow). 
    # Since we sorted descending, the curve is convex and decreasing.
    kneedle = KneeLocator(
        x, sorted_probs, 
        S=S, 
        curve="convex", 
        direction="decreasing"
    )
    
    knee_point = kneedle.knee
    
    # If no knee is found, return the top 5% as a fallback
    if knee_point is None:
        print("Warning: kneed found no elbow. Defaulting to top 5% features.")
        knee_point = max(1, int(len(sorted_probs) * 0.05))
        
    # Get the indices of the features up to the knee point
    top_indices = sorted_indices[:knee_point]
    
    # Map back to feature names
    key_features = [feature_names[i] for i in top_indices]
    
    return key_features
