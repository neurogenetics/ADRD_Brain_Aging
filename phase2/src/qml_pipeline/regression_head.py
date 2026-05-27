import torch
import torch.nn as nn

class AgeRegressionHead(nn.Module):
    """
    Unbounded continuous MLP regression head.
    Takes concatenated probabilities from the Walk Module and predicts biological age.
    """
    def __init__(self, input_dim: int, hidden_dim: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, 1) # No activation at the end for unbounded continuous prediction
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """
        Args:
            x: Concatenated, unpadded probabilities from the walk module (batch_size, total_genes)
        Returns:
            Predicted age (batch_size, 1)
        """
        # Ensure input is batched. The Walk module currently outputs a 1D tensor
        # for a single graph structure. In a real batched training loop, 
        # the graph topology is static, but we might batch over different patients.
        # If x is 1D (from a single forward pass), add batch dimension.
        if x.dim() == 1:
            x = x.unsqueeze(0)
            
        return self.net(x)

class CompletePipeline(nn.Module):
    """
    Combines the Walk Module and Regression Head into a single trainable model.
    """
    def __init__(self, walk_module: nn.Module):
        super().__init__()
        self.walk_module = walk_module
        
        # Calculate total input dimension for the MLP
        total_genes = sum(walk_module.original_node_counts.values())
        self.regression_head = AgeRegressionHead(input_dim=total_genes)
        
    def forward(self) -> torch.Tensor:
        # 1. Run the walk (classical or quantum)
        concat_probs = self.walk_module()
        
        # 2. Predict age
        age_pred = self.regression_head(concat_probs)
        return age_pred
