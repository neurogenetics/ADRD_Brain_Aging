import torch
import math
import numpy as np

def calculate_padding_size(num_nodes: int) -> int:
    """
    Calculates the number of ghost nodes required to pad a graph to 2^n nodes.
    This is a strict requirement for the Quantum Walk register.
    """
    if num_nodes <= 0:
        raise ValueError("Number of nodes must be greater than 0.")
    n_qubits = math.ceil(math.log2(num_nodes))
    target_nodes = 2 ** n_qubits
    return target_nodes - num_nodes

def compute_cosine_similarity(h_matrix: torch.Tensor) -> torch.Tensor:
    """
    Computes the dense cosine similarity matrix for a given cNMF H matrix.
    Args:
        h_matrix: Tensor of shape (num_factors, num_genes)
    Returns:
        Cosine similarity matrix of shape (num_genes, num_genes)
    """
    # Normalize along the factor dimension (L2 norm)
    h_norm = torch.nn.functional.normalize(h_matrix, p=2, dim=0)
    # Compute similarity (dot product of normalized vectors)
    sim_matrix = torch.mm(h_norm.T, h_norm)
    return sim_matrix

def pad_similarity_matrix(sim_matrix: torch.Tensor) -> torch.Tensor:
    """
    Pads a similarity matrix with ghost nodes (zeros) to reach 2^n dimension.
    """
    num_nodes = sim_matrix.shape[0]
    padding_size = calculate_padding_size(num_nodes)
    
    if padding_size == 0:
        return sim_matrix
        
    padded_matrix = torch.nn.functional.pad(
        sim_matrix, 
        (0, padding_size, 0, padding_size), 
        mode='constant', 
        value=0.0
    )
    return padded_matrix

def compute_graph_laplacian(sim_matrix: torch.Tensor, beta: torch.Tensor) -> torch.Tensor:
    """
    Computes the Graph Laplacian (L = D - A) using soft-thresholding.
    
    Args:
        sim_matrix: Padded cosine similarity matrix (2^n, 2^n)
        beta: Trainable parameter for soft-thresholding (must be positive)
    Returns:
        Graph Laplacian matrix
    """
    # Soft thresholding: A = |sim|^beta
    # Ensure no negative values in similarity (though cosine should be [-1, 1], 
    # cNMF H matrices are non-negative, so cosine sim is [0, 1])
    # Adding a small epsilon for numerical stability in backward pass if beta < 1 and sim=0
    adj_matrix = torch.pow(torch.abs(sim_matrix) + 1e-8, beta)
    
    # Ensure diagonal is zero (no self-loops in the adjacency matrix)
    # Using non-inplace multiplication with a mask to protect autograd
    diag_mask = 1.0 - torch.eye(adj_matrix.shape[0], device=adj_matrix.device, dtype=adj_matrix.dtype)
    adj_matrix = adj_matrix * diag_mask
    
    # Degree matrix D is a diagonal matrix where D_ii = sum(A_i)
    degree_vector = torch.sum(adj_matrix, dim=1)
    degree_matrix = torch.diag(degree_vector)
    
    # Graph Laplacian: L = D - A
    laplacian = degree_matrix - adj_matrix
    return laplacian

def sparsify_laplacian(laplacian: torch.Tensor, sparsity_threshold: float = 0.01) -> torch.Tensor:
    """
    Sets small values in the Laplacian to zero to create a sparse representation.
    """
    mask = torch.abs(laplacian) > sparsity_threshold
    # Keep the diagonal regardless of value
    diag_mask = torch.eye(laplacian.shape[0], device=laplacian.device, dtype=torch.bool)
    mask = mask | diag_mask
    return laplacian * mask
