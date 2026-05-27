import torch
import torch.nn as nn
import torch.optim as optim
from typing import List, Dict, Any
from .regression_head import CompletePipeline
import gc

def train_pipeline(
    model: CompletePipeline,
    target_age: float,
    num_epochs: int = 100,
    lr: float = 0.001,
    patience: int = 10
) -> Dict[str, Any]:
    """
    Trains the complete Walk + Regression pipeline.
    Uses HuberLoss to handle biological age outliers.
    """
    device = model.walk_module.device
    model.to(device)
    
    # Use gradient checkpointing inside the model to save memory
    # But since it's manual loops, we will just use torch.utils.checkpoint
    import torch.utils.checkpoint
    
    criterion = nn.HuberLoss()
    optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    
    target_tensor = torch.tensor([[target_age]], dtype=torch.float32, device=device)
    
    best_loss = float('inf')
    epochs_no_improve = 0
    history = []
    
    model.train()
    for epoch in range(num_epochs):
        torch.cuda.empty_cache()
        gc.collect()
        optimizer.zero_grad(set_to_none=True)
        
        # Use torch checkpoint to trade compute for memory
        def forward_wrapper():
            return model()
            
        predicted_age = torch.utils.checkpoint.checkpoint(forward_wrapper, use_reentrant=False)
        
        loss = criterion(predicted_age, target_tensor)
        loss.backward()
        
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()
        
        current_loss = loss.item()
        history.append(current_loss)
        
        if epoch % 10 == 0:
            print(f"Epoch {epoch}/{num_epochs} | Loss: {current_loss:.4f} | Pred: {predicted_age.item():.2f} | Target: {target_age}")
            
        if current_loss < best_loss:
            best_loss = current_loss
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= patience:
                print(f"Early stopping at epoch {epoch}")
                break
                
    return {
        "final_loss": history[-1],
        "history": history,
        "predicted_age": predicted_age.item()
    }
