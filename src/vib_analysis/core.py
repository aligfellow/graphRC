import numpy as np
from ase import Atoms
from itertools import combinations
import os
import logging
from typing import List, Dict, Tuple, Any, Union

from . import config

logger = logging.getLogger("vib_analysis")

def read_xyz_trajectory(file_path: str) -> List[Atoms]:
    """
    Read an XYZ trajectory file and return a list of ASE Atoms objects.
    
    Args:
        file_path: Path to XYZ trajectory file
        
    Returns:
        List of ASE Atoms objects, one per frame
        
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
            frame = Atoms(symbols=symbols, positions=coords)
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
    frame: Atoms, 
    bond_tolerance: float = config.BOND_TOLERANCE,
) -> Dict[str, List[Tuple[int, ...]]]:
    """
    Build internal coordinates (bonds, angles, dihedrals) using xyzgraph with hierarchical thresholds.
    
    Uses three separate graphs with progressively tighter thresholds:
    - Bond graph (flexible): captures forming/breaking bonds in TS
    - Tighter graph: more conservative connectivity for angles/dihedrals
    
    Args:
        frame: ASE Atoms object
        bond_tolerance: Multiplier for bond detection (most flexible)
        
    Returns:
        Dictionary with keys 'bonds', 'angles', 'dihedrals' containing lists of
        atom index tuples
    """
    from xyzgraph import build_graph
    
    # Convert ASE Atoms to xyzgraph format
    atoms_list = [(frame.get_chemical_symbols()[i], tuple(frame.positions[i])) 
                  for i in range(len(frame))]

    # xyzgraph base thresholds for flexible bond graph
    threshold_h_nm = 0.42 * bond_tolerance
    threshold_nm_nm = 0.55 * bond_tolerance 
    threshold_h_m = 0.50 * bond_tolerance
    threshold_m_l = 0.65  # Keep metal-ligand tight (no scaling)
    threshold_h_h = 0.40 * bond_tolerance 
    
    # 1. Build flexible graph for bonds
    G_bonds = build_graph(
        atoms=atoms_list,
        threshold=1.0,  # No additional global scaling
        threshold_h_nonmetal=threshold_h_nm,
        threshold_nonmetal_nonmetal=threshold_nm_nm,
        threshold_h_metal=threshold_h_m,
        threshold_metal_ligand=threshold_m_l,
        threshold_h_h=threshold_h_h,
        quick=True,  # Fast mode - we don't need bond orders
        method='cheminf',
        debug=False
    )
    
    # Extract bonds from flexible graph
    bonds = [(i, j) for i, j in G_bonds.edges()]
    logger.debug(f"Built bond graph with {len(bonds)} bonds (tolerance={bond_tolerance})")
    
    # 2. Build tighter graph for angles/dihedrals with default xyzgraph thresholds
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

def calculate_bond_length(frame: Atoms, i: int, j: int) -> float:
    """Calculate bond length between atoms i and j."""
    return round(float(frame.get_distance(i, j)), 3)


def calculate_angle(frame: Atoms, i: int, j: int, k: int) -> float:
    """Calculate angle between atoms i-j-k."""
    return round(float(frame.get_angle(i, j, k, mic=True)), 3)


def calculate_dihedral(frame: Atoms, i: int, j: int, k: int, l: int) -> float:
    """Calculate dihedral angle for atoms i-j-k-l."""
    return round(float(frame.get_dihedral(i, j, k, l, mic=True)), 3)


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
    frames: List[Atoms],
    ts_frame: Atoms,
    internal_coords: Dict[str, List[Tuple[int, ...]]],
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    bond_stability_threshold: float = config.BOND_STABILITY_THRESHOLD
) -> Tuple[Dict, Dict, Dict, Dict, Dict]:
    """
    Track changes in internal coordinates across trajectory.
    
    Identifies significant bond, angle, and dihedral changes between frames.
    Separates major changes from minor/dependent changes based on coupling
    to other structural changes.
    
    Args:
        frames: List of trajectory frames (typically 2 most diverse)
        ts_frame: Reference frame (typically transition state)
        internal_coords: Dictionary with 'bonds', 'angles', 'dihedrals' lists
        bond_threshold: Minimum bond change to report (Å)
        angle_threshold: Minimum angle change to report (degrees)
        dihedral_threshold: Minimum dihedral change to report (degrees)
        bond_stability_threshold: Threshold for filtering coupled angle/dihedral changes (Å)
        
    Returns:
        Tuple of (bond_changes, angle_changes, minor_angles, 
                  unique_dihedrals, dependent_dihedrals)
        Each is a dict mapping coordinate tuple to (max_change, initial_value)
    """
    # Identify significant bond changes
    bond_changes = {}
    for i, j in internal_coords['bonds']:
        distances = [calculate_bond_length(frame, i, j) for frame in frames]
        max_change = round(max(distances) - min(distances), 3)
        if abs(max_change) >= bond_threshold:
            initial_length = calculate_bond_length(ts_frame, i, j)
            bond_changes[(i, j)] = (max_change, initial_length)
    
    logger.debug(f"Found {len(bond_changes)} significant bond changes")
    
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
        
        # Skip if constituent bonds are not stable
        if not _bonds_are_stable(bonds_in_angle, bond_changes, bond_stability_threshold):
            continue
        
        # Calculate angle change
        angles = [calculate_angle(frame, i, j, k) for frame in frames]
        max_change = round(max(angles) - min(angles), 3)
        
        if abs(max_change) >= angle_threshold:
            initial_angle = calculate_angle(ts_frame, i, j, k)
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
        
        # Skip if constituent bonds are not stable
        if not _bonds_are_stable(bonds_in_dihedral, bond_changes, bond_stability_threshold):
            continue
        
        # Calculate dihedral change (adjust for periodicity)
        dihedrals = [calculate_dihedral(frame, i, j, k, l) for frame in frames]
        max_change = round(
            max([abs((d - dihedrals[0] + 180) % 360 - 180) for d in dihedrals]), 
            3
        )
        
        if max_change >= angle_threshold:
            dihedral_changes[(i, j, k, l)] = max_change
    
    # Group dihedrals by rotation axis and select representative
    masses = frames[0].get_masses()
    dihedral_groups = {}
    
    for (i, j, k, l), change in dihedral_changes.items():
        axis = tuple(sorted((j, k)))
        total_mass = masses[i] + masses[j] + masses[k] + masses[l]
        
        if axis not in dihedral_groups:
            dihedral_groups[axis] = []
        dihedral_groups[axis].append(((i, j, k, l), change, total_mass))
    
    # Select most significant dihedral per axis
    unique_dihedrals = {}
    dependent_dihedrals = {}
    
    for axis, dihedrals_list in dihedral_groups.items():
        # Sort by mass and change magnitude
        dihedrals_sorted = sorted(
            dihedrals_list, 
            key=lambda x: (x[2], x[1]), 
            reverse=True
        )
        dihedral, max_change, _ = dihedrals_sorted[0]
        
        if max_change >= dihedral_threshold:
            initial_dihedral = calculate_dihedral(ts_frame, *dihedral)
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
    
    return bond_changes, angle_changes, minor_angles, unique_dihedrals, dependent_dihedrals

def compute_rmsd(frame1, frame2):
    """Computes RMSD between two ASE frames."""
    diff = frame1.get_positions() - frame2.get_positions()
    return np.sqrt(np.mean(np.sum(diff**2, axis=1)))


def select_most_diverse_frames(frames, top_n=2):
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
    xyz_file_or_frames: Union[str, List[Atoms]],
    bond_tolerance: float = config.BOND_TOLERANCE,
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    bond_stability_threshold: float = config.BOND_STABILITY_THRESHOLD,
    ts_frame: int = config.DEFAULT_TS_FRAME,
) -> Dict[str, Any]:
    """
    Analyze vibrational displacements in trajectory to identify structural changes.
    
    This is the main analysis function that identifies bonds, angles, and dihedrals
    that change significantly across a trajectory. It selects the most diverse frames
    and tracks coordinate changes.
    
    Args:
        xyz_file_or_frames: Path to XYZ trajectory file or list of ASE Atoms objects
        bond_tolerance: Multiplier for xyzgraph bond detection thresholds
        bond_threshold: Minimum bond change to report (Å)
        angle_threshold: Minimum angle change to report (degrees)
        dihedral_threshold: Minimum dihedral change to report (degrees)
        bond_stability_threshold: Threshold for filtering coupled angle/dihedral changes (Å)
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
        TypeError: If xyz_file_or_frames is not str or list of Atoms
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
            "list of ASE Atoms objects."
        )

    internal_coords = build_internal_coordinates(
        frame=frames[ts_frame],
        bond_tolerance=bond_tolerance,
    )
    
    selected_indices = select_most_diverse_frames(frames)
    selected_frames = [frames[i] for i in selected_indices]
    
    logger.info(f"Selected frames {selected_indices} for analysis")

    bond_changes, angle_changes, minor_angles, unique_dihedrals, dependent_dihedrals = calculate_internal_changes(
        frames=selected_frames,
        ts_frame=frames[ts_frame],
        internal_coords=internal_coords,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        bond_stability_threshold=bond_stability_threshold,
    )

    first_frame = frames[0]
    symbols = first_frame.get_chemical_symbols()
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
    }
