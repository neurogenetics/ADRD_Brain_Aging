from typing import List

def jaccard_similarity(list1: List[str], list2: List[str]) -> float:
    """
    Computes the Jaccard similarity between two lists of features.
    Used to compare the overlap of modules extracted by Classical vs Quantum walks.
    """
    set1 = set(list1)
    set2 = set(list2)
    
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    
    if union == 0:
        return 0.0
        
    return intersection / union

def compare_walk_features(classical_features: List[str], quantum_features: List[str]) -> dict:
    """
    Generates a comparison report between the features extracted by the two walk methods.
    """
    j_sim = jaccard_similarity(classical_features, quantum_features)
    
    classical_only = list(set(classical_features) - set(quantum_features))
    quantum_only = list(set(quantum_features) - set(classical_features))
    intersection = list(set(classical_features).intersection(set(quantum_features)))
    
    return {
        "jaccard_similarity": j_sim,
        "num_classical_features": len(classical_features),
        "num_quantum_features": len(quantum_features),
        "num_shared": len(intersection),
        "classical_specific": classical_only,
        "quantum_specific": quantum_only
    }
