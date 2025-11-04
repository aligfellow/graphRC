"""
Configuration constants for vibrational analysis.

All default values can be overridden by passing parameters to analysis functions.
"""

# ============================================================================
# Coordinate Detection Tolerances
# ============================================================================
BOND_TOLERANCE = 1.4        # Multiplier for bond detection
RELAXED = False             # Use relaxed rules for xyzgraph

# ============================================================================
# Significance Thresholds
# ============================================================================
# Minimum changes to report as significant structural changes
BOND_THRESHOLD = 0.4        # Minimum bond change (Å)
ANGLE_THRESHOLD = 10.0      # Minimum angle change (degrees)
DIHEDRAL_THRESHOLD = 20.0   # Minimum dihedral change (degrees)

# Coupled proton transfer detection
COUPLED_PROTON_THRESHOLD = 0.15 # Low threshold for coupled H movements (Å)

# Internal filtering threshold
COUPLED_MOTION_FILTER  = 0.2  # For filtering coupled angle/dihedral changes

# ============================================================================
# Graph Analysis Parameters
# ============================================================================
DISTANCE_TOLERANCE = 0.2        # Bond formation/breaking tolerance (Å)
BOND_ORDER_EPSILON = 0.1        # Minimum change to report bond order change
IG_FLEXIBLE_DEFAULT = False     # Apply bond_tolerance to displaced graphs (with -ig)

# ============================================================================
# ASCII Visualization
# ============================================================================
ASCII_SCALE = 2.5               # Scaling factor for molecular rendering
ASCII_NEIGHBOR_SHELLS = 1       # Neighbor shells around reactive center
ASCII_INCLUDE_H = False         # Show hydrogen atoms by default

# ============================================================================
# Frame Selection
# ============================================================================
DEFAULT_TS_FRAME = 0            # Index of transition state frame
MAX_DISPLACEMENT_LEVEL = 4      # Maximum displacement amplitude level

# ============================================================================
# File Handling
# ============================================================================
SAVE_TRAJECTORY_DEFAULT = True
SAVE_DISPLACEMENT_DEFAULT = False
DEFAULT_DISPLACEMENT_LEVEL = 1
