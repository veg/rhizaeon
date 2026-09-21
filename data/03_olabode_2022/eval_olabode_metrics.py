"""
benchmarks/olabode_2022_benchmark/eval_olabode_metrics.py
========================================================
Standardized evaluation metrics following Olabode et al. (2022) PNAS:
1. Nucleotide-level partition misclassification error (%) via optimal Hungarian matching.
2. Breakpoint spatial error distance (nucleotides and codons).
3. Precision, recall, and detection accuracy under varying tolerances.
"""

import numpy as np
from scipy.optimize import linear_sum_assignment
from typing import List, Dict, Tuple, Optional


def build_partition_array(breakpoints: List[int], length: int) -> np.ndarray:
    """Constructs a 1D array of cluster/partition IDs across [0, length)."""
    bps = sorted([b for b in breakpoints if 0 < b < length])
    arr = np.zeros(length, dtype=int)
    prev = 0
    for cluster_id, b in enumerate(bps):
        arr[prev:b] = cluster_id
        prev = b
    arr[prev:length] = len(bps)
    return arr


def compute_partition_misclassification_error(
    true_bps: List[int],
    pred_bps: List[int],
    length: int
) -> float:
    """
    Computes nucleotide-level partition misclassification error (%)
    matching Olabode et al. (2022) Figure 1A.
    
    Constructs the contingency matrix between true and predicted partition labels
    and solves the optimal assignment (Hungarian algorithm) to find the maximum
    concordance.
    
    Returns:
        Error percentage in [0.0, 100.0]
    """
    if length <= 0:
        return 0.0
        
    true_labels = build_partition_array(true_bps, length)
    pred_labels = build_partition_array(pred_bps, length)
    
    n_true = len(true_bps) + 1
    n_pred = len(pred_bps) + 1
    
    # Compute contingency table
    # contingency[i, j] = number of positions where true_labels == i and pred_labels == j
    contingency = np.zeros((n_true, n_pred), dtype=np.int64)
    for x in range(length):
        contingency[true_labels[x], pred_labels[x]] += 1
        
    # Solve linear sum assignment on negative contingency matrix (maximize agreement)
    row_ind, col_ind = linear_sum_assignment(-contingency)
    matched_agreement = contingency[row_ind, col_ind].sum()
    
    misclassification_error = (1.0 - (matched_agreement / float(length))) * 100.0
    return float(misclassification_error)


def compute_breakpoint_spatial_metrics(
    true_bps: List[int],
    pred_bps: List[int],
    tolerances: List[int] = [30, 60, 150, 300]
) -> Dict[str, float]:
    """
    Computes spatial localization error (mean absolute distance) and
    hit rates under specified nucleotide tolerances.
    """
    if len(true_bps) == 0:
        return {
            "num_true": 0,
            "num_pred": len(pred_bps),
            "false_positives": len(pred_bps),
            "mean_error_nt": 0.0,
            "min_error_nt": 0.0,
            "recall_150nt": 1.0 if len(pred_bps) == 0 else 0.0,
            "precision_150nt": 1.0 if len(pred_bps) == 0 else 0.0
        }
        
    if len(pred_bps) == 0:
        res = {
            "num_true": len(true_bps),
            "num_pred": 0,
            "false_positives": 0,
            "mean_error_nt": float(np.nan),
            "median_error_nt": float(np.nan),
            "min_error_nt": float(np.nan),
        }
        for tol in tolerances:
            res[f"recall_{tol}nt"] = 0.0
            res[f"precision_{tol}nt"] = 0.0
        return res
        
    true_arr = np.array(sorted(true_bps))
    pred_arr = np.array(sorted(pred_bps))
    
    # For each true breakpoint, distance to closest predicted
    dist_true_to_pred = np.array([np.min(np.abs(pred_arr - tb)) for tb in true_arr])
    # For each predicted breakpoint, distance to closest true
    dist_pred_to_true = np.array([np.min(np.abs(true_arr - pb)) for pb in pred_arr])
    
    metrics = {
        "num_true": len(true_bps),
        "num_pred": len(pred_bps),
        "mean_error_nt": float(np.mean(dist_true_to_pred)),
        "median_error_nt": float(np.median(dist_true_to_pred)),
        "min_error_nt": float(np.min(dist_true_to_pred)),
        "max_error_nt": float(np.max(dist_true_to_pred))
    }
    
    for tol in tolerances:
        tp = np.sum(dist_true_to_pred <= tol)
        rec = tp / len(true_bps)
        fp = np.sum(dist_pred_to_true > tol)
        prec = (len(pred_bps) - fp) / len(pred_bps) if len(pred_bps) > 0 else 0.0
        metrics[f"recall_{tol}nt"] = float(rec)
        metrics[f"precision_{tol}nt"] = float(prec)
        
    return metrics


if __name__ == "__main__":
    # Self-test cases
    err0 = compute_partition_misclassification_error([4500], [4500], 9000)
    assert abs(err0 - 0.0) < 1e-6, f"Expected 0.0, got {err0}"
    
    err100 = compute_partition_misclassification_error([4500], [4600], 9000)
    assert abs(err100 - (100.0 / 9000 * 100)) < 1e-6, f"Expected 1.11%, got {err100}"
    
    err_miss = compute_partition_misclassification_error([4500], [], 9000)
    assert abs(err_miss - 50.0) < 1e-6, f"Expected 50.0%, got {err_miss}"
    
    print("All metric unit tests passed successfully!")
