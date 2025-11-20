import numpy as np
from itertools import combinations
import os
import logging
from typing import List, Dict, Tuple, Any, Union, Optional
from xyzgraph import DATA

from . import config
from .utils import calculate_distance, calculate_angle, calculate_dihedral

logger = logging.getLogger("vib_analysis")

def read_xyz_trajectory(file_path: str) -> List[Dict[str, Any]]:
    """
    Read an XYZ trajectory file and return a list of frame dicts.
    
    Args:
        file_path: Path to XYZ trajectory file
        
    Returns:
        List of frame dicts with 'symbols' and 'positions' keys, one per frame
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If only one geometry found (need at least 2 frames)
    """
    if not os.path.exists(file_path):
        logger.error(f"File {file_path} does not exist")
        raise FileNotFoundError(f"File {file_path} does not exist.") 
    
    frames = []
    with open(file_path, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            num_atoms = int(line.strip())
            _ = f.readline()  # Skip comment/title
            coords = []
            symbols = []
            for _ in range(num_atoms):
                parts = f.readline().split()
                symbols.append(parts[0])
                coords.append([float(x) for x in parts[1:]])
            frame = {
                'symbols': symbols,
                'positions': np.array(coords)
            }
            frames.append(frame)
    
    if len(frames) == 1:
        logger.error("Only one geometry found in trajectory file")
        raise ValueError(
            "Only one geometry found, make sure this is a trajectory file "
            "with at least 2 frames."
        )
    
    logger.info(f"Read {len(frames)} frames from {file_path}")
    return frames

def build_internal_coordinates(
    frame: Dict[str, Any],
    displaced_frames: Optional[List[Dict[str, Any]]] = None,
    independent_graphs: bool = False,
    ig_flexible: bool = config.IG_FLEXIBLE_DEFAULT,
    relaxed: bool = config.RELAXED,
    bond_tolerance: float = config.BOND_TOLERANCE,
) -> Dict[str, Any]:
    """
    Build internal coordinates (bonds, angles, dihedrals) using xyzgraph with hierarchical thresholds.
    Uses two separate graphs:
    - Bond graph (flexible): captures forming/breaking bonds in TS
    - Tighter graph: more conservative connectivity for angles/dihedrals

    When independent_graphs=True:
    - TS connectivity built with bond_tolerance (flexible)
    - Displaced connectivity built with xyzgraph defaults (strict) or bond_tolerance (if ig_flexible=True)
    - All connectivity merged via union
    - Bond coordinates derived from merged connectivity
    - Angles/dihedrals derived from tighter TS connectivity for consistency
    
    Args:
        frame: Frame dict with 'symbols' and 'positions' keys
        displaced_frames: Optional list of displaced frame dicts for connectivity augmentation
        independent_graphs: If True, augment TS connectivity with displaced connectivity
        ig_flexible: If True (with independent_graphs), apply bond_tolerance to displaced graphs
        relaxed: Use relaxed rules for xyzgraph (for complex ring systems)
        bond_tolerance: Multiplier for bond detection (most flexible)
        
    Returns:
        Dictionary with keys 'bonds', 'angles', 'dihedrals', 'connectivity'
    """
    from xyzgraph import build_graph
    
    # Convert frame dict to xyzgraph format
    symbols = frame['symbols']
    positions = frame['positions']
    atoms_list = [(symbols[i], tuple(positions[i])) 
                  for i in range(len(symbols))]

    # 1. Build flexible graph for bonds
    G_bonds = build_graph(
        atoms=atoms_list,
        threshold=bond_tolerance,  # No additional global scaling
        relaxed=relaxed,  # Use relaxed rules if specified
        quick=True,  # Fast mode - we don't need bond orders
        method='cheminf',
        debug=False
    )
    
    # Extract bonds from flexible graph
    ts_bonds = set(G_bonds.edges())
    logger.debug(f"Built TS bond graph with {len(ts_bonds)} bonds (tolerance={bond_tolerance})")
    
    # 2. If independent_graphs, augment with displaced connectivity
    if independent_graphs and displaced_frames:
        logger.debug(f"Augmenting TS connectivity with displaced graphs (flexible={ig_flexible})")
        
        augmented_bonds = set(ts_bonds)  # Start with TS
        
        disp_threshold = None
        # Determine thresholds for displaced graphs
        if ig_flexible:
            disp_threshold = bond_tolerance
            logger.debug("Using flexible thresholds for displaced graphs")
        else:
            disp_threshold = 1.0
            logger.debug("Using default thresholds for displaced graphs")
        
        # Build displaced graphs and collect bonds
        for idx, disp_frame in enumerate(displaced_frames):
            disp_atoms = [(disp_frame['symbols'][i], tuple(disp_frame['positions'][i])) 
                         for i in range(len(disp_frame['symbols']))]
            
            G_disp = build_graph(
                atoms=disp_atoms,
                threshold=disp_threshold,
                relaxed=relaxed,
                quick=True,
                method='cheminf',
                debug=False
            )
            
            disp_bonds = set(G_disp.edges())
            new_bonds = disp_bonds - ts_bonds
            if new_bonds:
                logger.debug(f"Displaced frame {idx}: found {len(new_bonds)} new bonds")
            augmented_bonds.update(disp_bonds)
        
        logger.debug(f"Augmented connectivity: {len(ts_bonds)} TS bonds + {len(augmented_bonds - ts_bonds)} new = {len(augmented_bonds)} total")
        bonds = list(augmented_bonds)
    else:
        # Standard behavior: TS only
        bonds = list(ts_bonds)
    
    # 3. Build tighter graph for angles/dihedrals with default xyzgraph thresholds
    # Always use TS geometry for consistency
    G_tight = build_graph(
        atoms=atoms_list,
        threshold=1.0,  # Use xyzgraph defaults
        quick=True,
        method='cheminf',
        debug=False
    )
    
    # Derive angles from tighter connectivity
    angles = []
    for j in G_tight.nodes():
        neighbors = list(G_tight.neighbors(j))
        if len(neighbors) >= 2:
            for i, k in combinations(neighbors, 2):
                angles.append((int(i), int(j), int(k)))
    
    # Derive dihedrals from tightest connectivity
    dihedrals = []
    for b, c in G_tight.edges():
        a_neighbors = set(G_tight.neighbors(b)) - {c}
        d_neighbors = set(G_tight.neighbors(c)) - {b}

        for a in a_neighbors:
            for d in d_neighbors:
                if a != d:
                    dihedrals.append((int(a), int(b), int(c), int(d)))
    
    # Build connectivity dictionary from tight graph for angle/dihedral consistency
    # Use G_tight since dihedrals (used for inversion detection) come from this graph
    connectivity = {}
    for i, j in G_tight.edges():
        connectivity.setdefault(i, set()).add(j)
        connectivity.setdefault(j, set()).add(i)
    
    return {
        'bonds': bonds, 
        'angles': angles, 
        'dihedrals': dihedrals,
        'connectivity': connectivity
    }

def _has_significant_bond_change(
    bond: Tuple[int, int], 
    bond_changes: Dict[Tuple[int, int], Tuple[float, float]], 
    threshold: float
) -> bool:
    """
    Check if a bond has a significant change above threshold.
    
    Args:
        bond: Tuple of atom indices
        bond_changes: Dictionary of bond changes
        threshold: Minimum change threshold
        
    Returns:
        True if bond change exceeds threshold
    """
    sorted_bond = tuple(sorted(bond))
    change_data = bond_changes.get(sorted_bond, (0.0, 0.0))
    return change_data[0] >= threshold


def _bonds_are_stable(
    bonds: List[Tuple[int, int]], 
    bond_changes: Dict[Tuple[int, int], Tuple[float, float]], 
    threshold: float
) -> bool:
    """
    Check if all bonds in list are below stability threshold.
    
    Args:
        bonds: List of bond tuples
        bond_changes: Dictionary of bond changes
        threshold: Stability threshold
        
    Returns:
        True if all bonds are stable (below threshold)
    """
    return all(
        bond_changes.get(tuple(sorted(bond)), (0.0, 0.0))[0] < threshold 
        for bond in bonds
    )


def calculate_internal_changes(
    frames: List[Dict[str, Any]],
    ts_frame: Dict[str, Any],
    internal_coords: Dict[str, Any],
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    coupled_motion_filter: float = config.COUPLED_MOTION_FILTER,
    coupled_proton_threshold: Union[float, bool] = config.COUPLED_PROTON_THRESHOLD
) -> Tuple[Dict, Dict, Dict, Dict, Dict, set]:
    """
    Track changes in internal coordinates across trajectory.
    
    Identifies significant bond, angle, and dihedral changes between frames.
    Separates major changes from minor/dependent changes based on coupling
    to other structural changes.
    
    Args:
        frames: List of frame dicts (typically 2 most diverse)
        ts_frame: Reference frame dict (typically transition state)
        internal_coords: Dictionary with 'bonds', 'angles', 'dihedrals' lists
        bond_threshold: Minimum bond change to report (Å)
        angle_threshold: Minimum angle change to report (degrees)
        dihedral_threshold: Minimum dihedral change to report (degrees)
        coupled_motion_filter: Threshold for filtering coupled angle/dihedral changes (Å)
        
    Returns:
        Tuple of (bond_changes, angle_changes, minor_angles, 
                  unique_dihedrals, dependent_dihedrals, coupled_proton_bonds)
        Each dict maps coordinate tuple to (max_change, initial_value)
        coupled_proton_bonds is a set of bonds detected via coupled threshold
    """
    # Get symbols from ts_frame for element identification
    symbols = ts_frame['symbols']
    
    # Identify significant bond changes
    bond_changes = {}
    for i, j in internal_coords['bonds']:
        distances = [calculate_distance(frame['positions'], i, j) for frame in frames]
        max_change = round(max(distances) - min(distances), 3)
        if abs(max_change) >= bond_threshold:
            initial_length = calculate_distance(ts_frame['positions'], i, j)
            bond_changes[(i, j)] = (max_change, initial_length)
    
    logger.debug(f"Found {len(bond_changes)} significant bond changes")
    
    # Coupled proton transfer detection (for any H involved in detected bond changes)
    coupled_proton_bonds = set()  # Track which bonds were detected via coupled threshold
    if coupled_proton_threshold is not False:
        # Find ALL H atoms involved in any detected bond change
        h_atoms_in_changes = set()
        for (i, j) in bond_changes:
            if symbols[i] == 'H':
                h_atoms_in_changes.add(i)
            if symbols[j] == 'H':
                h_atoms_in_changes.add(j)
        
        if h_atoms_in_changes:
            logger.debug(f"Searching for coupled proton transfers for {len(h_atoms_in_changes)} H atoms")
        
        # Search for additional bonds involving these H atoms with reduced threshold
        for h_idx in h_atoms_in_changes:
            for i, j in internal_coords['bonds']:
                if h_idx not in (i, j):
                    continue  # Skip bonds not involving this H
                
                sorted_bond = tuple(sorted((i, j)))
                if sorted_bond in bond_changes:
                    continue  # Already detected
                
                # Check with reduced threshold
                distances = [calculate_distance(frame['positions'], i, j) for frame in frames]
                max_change = round(max(distances) - min(distances), 3)
                if abs(max_change) >= coupled_proton_threshold:
                    initial_length = calculate_distance(ts_frame['positions'], i, j)
                    bond_changes[sorted_bond] = (max_change, initial_length)  
                    coupled_proton_bonds.add(sorted_bond)  # Track separately
                    logger.debug(f"Coupled proton detection: bond {sorted_bond} (Δ={max_change} Å)")
                    # break     # commented out to allow other checks of connectivity
    
    # Track atoms involved in bond changes
    changed_atoms = set()
    for bond in bond_changes:
        changed_atoms.update(bond)
    
    # Process angle changes
    angle_changes = {}
    minor_angles = {}
    
    for i, j, k in internal_coords['angles']:
        bonds_in_angle = [tuple(sorted((i, j))), tuple(sorted((j, k)))]
        
        # Skip if any constituent bond has significant change
        if any(bond in bond_changes for bond in bonds_in_angle):
            continue
        
        # Skip if bonds are not stable (some motion but not enough for reporting a bond change)
        if not _bonds_are_stable(bonds_in_angle, bond_changes, coupled_motion_filter):
            continue
        
        # Calculate angle change
        angles = [calculate_angle(frame['positions'], i, j, k) for frame in frames]
        max_change = round(max(angles) - min(angles), 3)
        
        if abs(max_change) >= angle_threshold:
            initial_angle = calculate_angle(ts_frame['positions'], i, j, k)
            angle_atoms = set((i, j, k))
            
            # Classify as minor if involves atoms from bond changes
            if angle_atoms.intersection(changed_atoms):
                minor_angles[(i, j, k)] = (max_change, initial_angle)
            else:
                angle_changes[(i, j, k)] = (max_change, initial_angle)
    
    logger.debug(
        f"Found {len(angle_changes)} significant angles, "
        f"{len(minor_angles)} minor angles"
    )
    
    # Process dihedral changes
    dihedral_changes = {}
    
    for i, j, k, l in internal_coords['dihedrals']:
        bonds_in_dihedral = [(i, j), (j, k), (k, l)]
        
        # Skip if any bond in dihedral has significant change
        if any(set(bond).issubset({i, j, k, l}) for bond in bond_changes):
            continue
        
        # Skip if bonds are not stable (some motion but not enough for reporting a bond change)
        if not _bonds_are_stable(bonds_in_dihedral, bond_changes, coupled_motion_filter):
            continue
        
        # Calculate dihedral change (adjust for periodicity)
        dihedrals = [calculate_dihedral(frame['positions'], i, j, k, l) for frame in frames]
        max_change = round(
            max([abs((d - dihedrals[0] + 180) % 360 - 180) for d in dihedrals]), 
            3
        )
        
        if max_change >= angle_threshold:
            dihedral_changes[(i, j, k, l)] = max_change
    
    # Group dihedrals by rotation axis and select representative
    # Use atomic numbers as proxy for mass (heavier atoms have higher atomic numbers)
    symbols = frames[0]['symbols']
    atomic_numbers = [DATA.s2n[sym] for sym in symbols]
    dihedral_groups = {}
    
    for (i, j, k, l), change in dihedral_changes.items():
        axis = tuple(sorted((j, k)))
        total_atomic_number = atomic_numbers[i] + atomic_numbers[j] + atomic_numbers[k] + atomic_numbers[l]
        
        if axis not in dihedral_groups:
            dihedral_groups[axis] = []
        dihedral_groups[axis].append(((i, j, k, l), change, total_atomic_number))
    
    # Select most significant dihedral per axis
    unique_dihedrals = {}
    dependent_dihedrals = {}
    
    for axis, dihedrals_list in dihedral_groups.items():
        # Sort by atomic number sum and change magnitude
        dihedrals_sorted = sorted(
            dihedrals_list, 
            key=lambda x: (x[2], x[1]), 
            reverse=True
        )
        dihedral, max_change, _ = dihedrals_sorted[0]
        
        if max_change >= dihedral_threshold:
            initial_dihedral = calculate_dihedral(ts_frame['positions'], *dihedral)
            dihedral_atoms = set(dihedral)
            
            # Classify as dependent if involves atoms from bond changes
            if dihedral_atoms.intersection(changed_atoms):
                dependent_dihedrals[dihedral] = (max_change, initial_dihedral)
            else:
                unique_dihedrals[dihedral] = (max_change, initial_dihedral)
    
    logger.debug(
        f"Found {len(unique_dihedrals)} significant dihedrals, "
        f"{len(dependent_dihedrals)} dependent dihedrals"
    )
    
    if coupled_proton_bonds:
        logger.debug(f"Found {len(coupled_proton_bonds)} bonds via coupled proton threshold")
    
    return bond_changes, angle_changes, minor_angles, unique_dihedrals, dependent_dihedrals, coupled_proton_bonds

def compute_rmsd(frame1: Dict[str, Any], frame2: Dict[str, Any]) -> float:
    """Computes RMSD between two frame dicts."""
    diff = frame1['positions'] - frame2['positions']
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


def select_most_diverse_frames(frames: List[Dict[str, Any]], top_n: int = 2) -> List[int]:
    """Select frames with largest RMSD from the TS frame (frame 0)."""
    # create an RMSD matrix between all frames and select the highest pair
    rmsd_matrix = np.zeros((len(frames), len(frames)))
    highest_rmsd = 0.0
    indices = []
    for i in range(len(frames)):
        for j in range(i + 1, len(frames)):
            rmsd_value = compute_rmsd(frames[i], frames[j])
            rmsd_matrix[i][j] = rmsd_value
            rmsd_matrix[j][i] = rmsd_value
    
    # get the largest RMSD value from the matrix
    for i in range(len(frames)):
        for j in range(i + 1, len(frames)):
            if rmsd_matrix[i][j] > highest_rmsd:
                highest_rmsd = rmsd_matrix[i][j]
                indices = [i, j]

    selected_indices = indices
    return selected_indices

def analyze_internal_displacements(
    xyz_file_or_frames: Union[str, List[Dict[str, Any]]],
    relaxed: bool = config.RELAXED,
    bond_tolerance: float = config.BOND_TOLERANCE,
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    coupled_motion_filter: float = config.COUPLED_MOTION_FILTER,
    coupled_proton_threshold: Union[float, bool] = config.COUPLED_PROTON_THRESHOLD,
    ts_frame: int = config.DEFAULT_TS_FRAME,
    independent_graphs: bool = False,
    ig_flexible: bool = config.IG_FLEXIBLE_DEFAULT,
) -> Dict[str, Any]:
    """
    Analyze vibrational displacements in trajectory to identify structural changes.
    
    This is the main analysis function that identifies bonds, angles, and dihedrals
    that change significantly across a trajectory. It selects the most diverse frames
    and tracks coordinate changes.
    
    Args:
        xyz_file_or_frames: Path to XYZ trajectory file or list of frame dicts
        bond_tolerance: Multiplier for xyzgraph bond detection thresholds
        bond_threshold: Minimum bond change to report (Å)
        angle_threshold: Minimum angle change to report (degrees)
        dihedral_threshold: Minimum dihedral change to report (degrees)
        coupled_motion_filter: Threshold for filtering coupled angle/dihedral changes (Å)
        ts_frame: Index of transition state frame to use as reference
        
    Returns:
        Dictionary containing:
            - bond_changes: Dict mapping bond tuples to (change, initial_value)
            - angle_changes: Dict mapping angle tuples to (change, initial_value)
            - minor_angle_changes: Dict of angles coupled to bond changes
            - dihedral_changes: Dict mapping dihedral tuples to (change, initial_value)
            - minor_dihedral_changes: Dict of dihedrals coupled to bond changes
            - frame_indices: List of selected frame indices
            - atom_index_map: Dict mapping atom indices to symbols
            
    Raises:
        TypeError: If xyz_file_or_frames is not str or list of frame dicts
        FileNotFoundError: If file path doesn't exist
        ValueError: If trajectory has less than 2 frames
    """
    # Handle both file path and in-memory frames
    if isinstance(xyz_file_or_frames, str):
        frames = read_xyz_trajectory(xyz_file_or_frames)
    elif isinstance(xyz_file_or_frames, list):
        frames = xyz_file_or_frames
    else:
        raise TypeError(
            "xyz_file_or_frames must be a file path (str) or "
            "list of frame dicts."
        )

    # Select diverse frames FIRST (needed for independent_graphs mode)
    selected_indices = select_most_diverse_frames(frames)
    selected_frames = [frames[i] for i in selected_indices]
    
    logger.info(f"Using TS frame {ts_frame}, selected frames {selected_indices} for analysis")
    
    # Build internal coordinates with optional augmentation
    internal_coords = build_internal_coordinates(
        frame=frames[ts_frame],
        displaced_frames=selected_frames if independent_graphs else None,
        independent_graphs=independent_graphs,
        ig_flexible=ig_flexible,
        relaxed=relaxed,
        bond_tolerance=bond_tolerance,
    )

    bond_changes, angle_changes, minor_angles, unique_dihedrals, dependent_dihedrals, coupled_proton_bonds = calculate_internal_changes(
        frames=selected_frames,
        ts_frame=frames[ts_frame],
        internal_coords=internal_coords,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        coupled_motion_filter=coupled_motion_filter,
        coupled_proton_threshold=coupled_proton_threshold,
    )

    first_frame = frames[0]
    symbols = first_frame['symbols']
    atom_index_map = {i: s for i, s in enumerate(symbols)}

    return {
        "bond_changes": bond_changes,
        "angle_changes": angle_changes,
        "minor_angle_changes": minor_angles,
        "dihedral_changes": unique_dihedrals,
        "minor_dihedral_changes": dependent_dihedrals,
        "frame_indices": selected_indices,
        "atom_index_map": atom_index_map,
        "connectivity": internal_coords['connectivity'],
        "coupled_proton_bonds": coupled_proton_bonds,
    }
