"""
posebench-fast: Fast docking evaluation metrics

Provides:
- Symmetry-corrected RMSD computation
- Fast PoseBusters filters (without full energy evaluation)
- Docking metrics (success rates, averages)
"""

from posebench_fast.metrics.rmsd import (
    compute_all_isomorphisms,
    get_symmetry_rmsd_with_isomorphisms,
    get_symmetry_rmsd,
    symmrmsd,
    TimeoutException,
    time_limit,
)

from posebench_fast.filters.fast_filters import (
    calc_posebusters,
    check_intermolecular_distance,
    check_volume_overlap,
    check_geometry,
)

from posebench_fast.metrics.aggregation import (
    get_simple_metrics_df,
    get_final_results_for_df,
    filter_results_by_posebusters,
    filter_results_by_fast,
    get_best_results_by_score,
)

__version__ = "0.1.0"

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
