"""
Módulo de generación de trayectorias para robot 2R.
"""
from .trefoil import (
    generate_trefoil_raw,
    parametrize_constant_speed,
    find_tangential_entry,
    generate_bezier_approach,
    filter_savgol,
)
from .builder import build_full_trajectory
