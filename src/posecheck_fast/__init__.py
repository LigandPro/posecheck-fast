"""
posecheck-fast: Fast docking evaluation metrics

Provides:
- Symmetry-corrected RMSD computation
- Fast PoseBusters filters (without full energy evaluation)
- Docking metrics (success rates, averages)
"""

from posecheck_fast.filters.fast_filters import (
    calc_posebusters,
    check_geometry,
    check_intermolecular_distance,
    check_volume_overlap,
)
from posecheck_fast.metrics.aggregation import (
    filter_results_by_fast,
    filter_results_by_posebusters,
    get_best_results_by_score,
    get_final_results_for_df,
    get_simple_metrics_df,
)
from posecheck_fast.metrics.rmsd import (
    TimeoutException,
    compute_all_isomorphisms,
    get_symmetry_rmsd,
    get_symmetry_rmsd_with_isomorphisms,
    symmrmsd,
    time_limit,
)

__version__ = "0.1.12"
__all__ = [
    # RMSD
    "compute_all_isomorphisms",
    "get_symmetry_rmsd_with_isomorphisms",
    "get_symmetry_rmsd",
    "symmrmsd",
    "TimeoutException",
    "time_limit",
    # Filters
    "calc_posebusters",
    "check_intermolecular_distance",
    "check_volume_overlap",
    "check_geometry",
    # Metrics
    "get_simple_metrics_df",
    "get_final_results_for_df",
    "filter_results_by_posebusters",
    "filter_results_by_fast",
    "get_best_results_by_score",
]
