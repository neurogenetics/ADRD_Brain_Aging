import torch
import torch.nn as nn
import math
from typing import Dict, Optional, Tuple
from .graph_utils import compute_cosine_similarity, pad_similarity_matrix, compute_graph_laplacian, sparsify_laplacian

def split_mat_vec_mult(A: torch.Tensor, v: torch.Tensor, chunks: int = 4) -> torch.Tensor:
    """
    Performs chunked matrix-vector multiplication to save peak memory
    during the autograd backward pass.
    """
    chunk_size = A.shape[0] // chunks
    result = []
    for i in range(chunks):
        start = i * chunk_size
        end = start + chunk_size if i < chunks - 1 else A.shape[0]
        # Multiply a chunk of rows by the vector
        res_chunk = torch.matmul(A[start:end, :], v)
        result.append(res_chunk)
    return torch.cat(result, dim=0)

def taylor_matrix_exp_vector(A: torch.Tensor, v: torch.Tensor, order: int = 5) -> torch.Tensor:
    """
    Computes exp(A) * v using a Taylor series expansion.
    """
    result = v.clone()
    term = v.clone()
    
    # We clear the cache to ensure we have the max contiguous block available for matmuls
    torch.cuda.empty_cache()
    
    for k in range(1, order + 1):
        # term = torch.matmul(A, term) / k
        term = split_mat_vec_mult(A, term, chunks=4) / k
        result = result + term
        
        # Early stopping if the term becomes negligibly small
        if torch.max(torch.abs(term)) < 1e-4:
            break
            
    return result

class BaseWalkModule(nn.Module):
    """
    Base module for Continuous-Time Walks.
    Handles graph construction and VRAM-safe sequential processing.
    """
    def __init__(self, h_matrices: Dict[str, torch.Tensor], walk_type: str = "classical", device: str = "cuda"):
        super().__init__()
        self.walk_type = walk_type
        self.device = torch.device(device if torch.cuda.is_available() else "cpu")
        
        # Store original node counts before padding to strip padding later
        self.original_node_counts = {name: h.shape[1] for name, h in h_matrices.items()}
        
        # Trainable parameters
        self.beta = nn.Parameter(torch.tensor(2.0, dtype=torch.float32))
        self.time_step = nn.Parameter(torch.tensor(1.0, dtype=torch.float32))
        
        # Precompute and register padded similarity matrices to save computation
        self.sim_matrices = nn.ParameterDict()
        for name, h_mat in h_matrices.items():
            h_mat_dev = h_mat.to(self.device)
            sim = compute_cosine_similarity(h_mat_dev)
            padded_sim = pad_similarity_matrix(sim)
            # Use fp16 to save memory
            self.register_buffer(f"sim_{name}", padded_sim.half())

        self.names = list(h_matrices.keys())

    def forward(self) -> torch.Tensor:
        """
        Executes the walk sequentially across all inputs to protect VRAM.
        Returns a concatenated vector of unpadded probabilities.
        """
        # Ensure beta and time_step remain positive
        beta_pos = torch.nn.functional.softplus(self.beta)
        t_pos = torch.nn.functional.softplus(self.time_step)
        
        unpadded_probs = []
        
        for name in self.names:
            sim_matrix = getattr(self, f"sim_{name}").float() # cast back to float32 for autograd
            laplacian = compute_graph_laplacian(sim_matrix, beta_pos)
            
            # Optional: Sparsify to speed up matmuls if the graph is mostly disconnected after soft thresholding
            laplacian = sparsify_laplacian(laplacian, sparsity_threshold=1e-4)
            
            # Prevent exploding gradients by clipping Laplacian values
            laplacian = torch.clamp(laplacian, min=-10.0, max=10.0)
            
            num_nodes = laplacian.shape[0]
            
            if self.walk_type == "classical":
                # CTRW: p(t) = exp(-L * t) * p(0)
                initial_state = torch.ones(num_nodes, dtype=torch.float32, device=self.device) / num_nodes
                
                A = -laplacian * t_pos
                # order 5 is enough to get gradients flowing, reducing memory pressure 
                # (10 copies of the gradient graph is too much for 16GB VRAM on 16384^2)
                probs = taylor_matrix_exp_vector(A, initial_state, order=5)
                
            elif self.walk_type == "quantum":
                # CTQW: |psi(t)> = exp(-i * L * t) |psi(0)>
                laplacian_c = laplacian.to(torch.complex64)
                t_pos_c = t_pos.to(torch.complex64)
                
                initial_state = torch.ones(num_nodes, dtype=torch.complex64, device=self.device) / math.sqrt(num_nodes)
                
                A = -1j * laplacian_c * t_pos_c
                # Reduced order
                evolved_state = taylor_matrix_exp_vector(A, initial_state, order=5)
                
                probs = torch.abs(evolved_state) ** 2
            else:
                raise ValueError(f"Unknown walk type: {self.walk_type}")
                
            orig_count = self.original_node_counts[name]
            bio_probs = probs[:orig_count]
            
            # Add small epsilon to prevent NaN during log/backward
            bio_probs = bio_probs + 1e-8
            # Re-normalize just in case Taylor drifted
            bio_probs = bio_probs / torch.sum(bio_probs)
            
            unpadded_probs.append(bio_probs)
            
        return torch.cat(unpadded_probs, dim=0)

class ClassicalWalk(BaseWalkModule):
    def __init__(self, h_matrices: Dict[str, torch.Tensor], device: str = "cuda"):
        super().__init__(h_matrices, walk_type="classical", device=device)

class QuantumWalk(BaseWalkModule):
    def __init__(self, h_matrices: Dict[str, torch.Tensor], device: str = "cuda"):
        super().__init__(h_matrices, walk_type="quantum", device=device)
