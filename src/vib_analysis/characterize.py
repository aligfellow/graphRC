"""
Motion characterization for vibrational modes using simple, robust logic.

Uses topology-based analysis:
- Inversions: Identified by "hub" pattern where one atom dominates dihedral centers
- Rotations: Classified by neighbor counting (Me, tBu, etc.)

No complex geometry calculations - just counting and simple logic.
"""

import numpy as np
from typing import List, Dict, Tuple, Any, Optional
from ase import Atoms
import logging

logger = logging.getLogger("vib_analysis")


# ============================================================================
# DISPLACEMENT UTILITIES
# ============================================================================

def calculate_atom_displacement(
    atom_idx: int,
    frame_ts: Atoms,
    frames_displaced: List[Atoms]
) -> float:
    """
    Calculate average displacement magnitude for an atom.
    
    Args:
        atom_idx: Atom index
        frame_ts: TS frame
        frames_displaced: Displaced frames
    
    Returns:
        Average displacement magnitude in Angstroms
    """
    pos_ts = frame_ts.get_positions()[atom_idx]
    displacements = []
    
    for frame in frames_displaced:
        pos = frame.get_positions()[atom_idx]
        disp = np.linalg.norm(pos - pos_ts)
        displacements.append(disp)
    
    return np.mean(displacements) if displacements else 0.0


# ============================================================================
# HUB DETECTION FOR INVERSIONS
# ============================================================================

def detect_inversion_hub(
    dihedral_changes: Dict[Tuple[int, int, int, int], Tuple[float, float]],
    symbols: List[str],
    min_hub_fraction: float = 0.5
) -> Optional[Tuple[int, str, float]]:
    """
    Detect if dihedrals share a common central "hub" atom, indicating inversion.
    
    Inversions create a hub pattern where one atom appears as the central atom
    (positions j or k) in multiple dihedrals. For example, N inversion shows N
    in the center of all affected dihedrals.
    
    Args:
        dihedral_changes: Dict of (i,j,k,l) -> (max_change, initial_value)
        symbols: List of atomic symbols
        min_hub_fraction: Minimum fraction of dihedrals hub must appear in
    
    Returns:
        (hub_atom_index, element_symbol, hub_fraction) if inversion detected, else None
    """
    if len(dihedral_changes) < 3:
        # Need at least 3 dihedrals for reliable hub detection
        # With only 2 dihedrals, can't distinguish inversion from rotation
        return None
    
    # Count occurrences of each atom in central positions (j and k)
    hub_counts = {}
    for (i, j, k, l) in dihedral_changes.keys():
        hub_counts[j] = hub_counts.get(j, 0) + 1
        hub_counts[k] = hub_counts.get(k, 0) + 1
    
    if not hub_counts:
        return None
    
    # Find atom with highest count
    hub_atom = max(hub_counts, key=hub_counts.get)
    hub_count = hub_counts[hub_atom]
    hub_fraction = hub_count / len(dihedral_changes)
    
    logger.debug(f"Hub analysis: atom {hub_atom} ({symbols[hub_atom]}) appears in {hub_count}/{len(dihedral_changes)} dihedrals ({hub_fraction:.1%})")
    
    # Check if hub dominates (appears in at least min_hub_fraction of dihedrals)
    if hub_fraction >= min_hub_fraction:
        # Check if it's an invertible atom
        invertible_elements = {'N', 'P', 'S', 'As', 'Sb'}
        element = symbols[hub_atom]
        
        if element in invertible_elements:
            logger.info(
                f"Inversion detected: {element} at position {hub_atom} is hub for {hub_fraction:.1%} of dihedrals"
            )
            return (hub_atom, element, hub_fraction)
        else:
            logger.debug(f"Hub detected but not invertible element: {element}")
    
    return None


def identify_moving_group(
    hub_atom: int,
    frame_ts: Atoms,
    frames_displaced: List[Atoms],
    connectivity: Dict[int, set]
) -> Tuple[int, float, str]:
    """
    Identify which group on the hub atom moves most during inversion.
    
    Args:
        hub_atom: Index of inverting atom
        frame_ts: TS frame
        frames_displaced: Displaced frames
        connectivity: Connectivity dictionary from build_internal_coordinates
    
    Returns:
        (neighbor_idx, displacement, description) of most mobile group
    """
    # Get neighbors from connectivity
    neighbors = list(connectivity.get(hub_atom, set()))
    symbols = frame_ts.get_chemical_symbols()
    
    if len(neighbors) == 0:
        return (hub_atom, 0.0, "unknown")
    
    # Calculate displacement for each neighbor
    neighbor_displacements = []
    for neighbor in neighbors:
        disp = calculate_atom_displacement(neighbor, frame_ts, frames_displaced)
        neighbor_displacements.append((neighbor, disp, symbols[neighbor]))
    
    # Sort by displacement
    neighbor_displacements.sort(key=lambda x: x[1], reverse=True)
    
    most_mobile = neighbor_displacements[0]
    neighbor_idx, max_disp, neighbor_symbol = most_mobile
    
    # Try to identify the group type
    if neighbor_symbol == 'C':
        # Check if it's part of a methyl (C with 3H neighbors)
        c_neighbors = connectivity.get(neighbor_idx, set())
        h_count = sum(1 for n in c_neighbors if symbols[n] == 'H')
        if h_count == 3:
            description = "methyl group"
        else:
            description = f"{neighbor_symbol} group"
    elif neighbor_symbol == 'H':
        description = "hydrogen"
    else:
        description = f"{neighbor_symbol} group"
    
    logger.debug(
        f"Most mobile group: {description} at position {neighbor_idx} "
        f"(displacement: {max_disp:.3f} Å)"
    )
    
    return (neighbor_idx, max_disp, description)


# ============================================================================
# ROTATION CLASSIFICATION
# ============================================================================

def find_aromatic_rings(frame: Atoms, connectivity: Dict[int, set]) -> List[List[int]]:
    """
    Find aromatic 6-membered carbon rings using simple heuristic.
    
    Criteria:
    - 6-membered ring
    - All carbons
    - Each carbon has 3-5 neighbors (sp2 + substituents)
    
    Args:
        frame: ASE Atoms object
        connectivity: Connectivity dictionary from build_internal_coordinates
    
    Returns:
        List of rings, where each ring is a list of atom indices
    """
    symbols = frame.get_chemical_symbols()
    aromatic_rings = []
    
    # Find all 6-membered rings using BFS
    from collections import deque
    
    visited_rings = set()
    
    for start_atom in range(len(frame)):
        if symbols[start_atom] != 'C':
            continue
        
        # BFS to find 6-membered rings starting from this atom
        queue = deque([(start_atom, [start_atom], {start_atom})])
        
        while queue:
            current, path, visited = queue.popleft()
            
            if len(path) > 6:
                continue
            
            neighbors = list(connectivity.get(current, set()))
            for next_atom in neighbors:
                # Check if we've completed a 6-membered ring
                if next_atom == start_atom and len(path) == 6:
                    # Found 6-membered ring - check if aromatic
                    ring_tuple = tuple(sorted(path))
                    if ring_tuple not in visited_rings:
                        # Check: all carbons with 3-5 neighbors (allowing substituents)
                        is_aromatic = True
                        for atom in path:
                            if symbols[atom] != 'C':
                                is_aromatic = False
                                break
                            atom_neighbors = connectivity.get(atom, set())
                            # Allow 3-5 neighbors (ring + 0-3 substituents)
                            if len(atom_neighbors) < 3 or len(atom_neighbors) > 5:
                                is_aromatic = False
                                break
                        
                        if is_aromatic:
                            aromatic_rings.append(path)
                            visited_rings.add(ring_tuple)
                
                # Continue search
                if next_atom not in visited and len(path) < 6:
                    new_visited = visited.copy()
                    new_visited.add(next_atom)
                    queue.append((next_atom, path + [next_atom], new_visited))
    
    return aromatic_rings


def classify_rotation_type(
    dihedral: Tuple[int, int, int, int],
    frame_ts: Atoms,
    dihedral_change: float,
    connectivity: Dict[int, set]
) -> Dict[str, Any]:
    """
    Classify type of rotation based on chemical environment.
    
    Conservative logic - only classify when confident:
    - Check for aromatic ring membership
    - Only classify symmetric tops (Me) when not in rings
    - Default to generic descriptions otherwise
    
    Args:
        dihedral: (i, j, k, l) atom indices
        frame_ts: TS frame
        dihedral_change: Rotation magnitude (degrees)
        connectivity: Connectivity dictionary from build_internal_coordinates
    
    Returns:
        Dict with rotation characterization
    """
    i, j, k, l = dihedral
    symbols = frame_ts.get_chemical_symbols()
    
    # Get neighbors from connectivity
    neighbors_j = list(connectivity.get(j, set()))
    neighbors_k = list(connectivity.get(k, set()))
    
    # Find aromatic rings
    aromatic_rings = find_aromatic_rings(frame_ts, connectivity)
    
    # Check if any dihedral atoms are in aromatic rings
    aromatic_atoms = []
    for ring in aromatic_rings:
        ring_set = set(ring)
        if j in ring_set:
            aromatic_atoms.append(j)
        if k in ring_set:
            aromatic_atoms.append(k)
    
    # If aromatic atoms detected, report aromatic rotation
    if aromatic_atoms:
        aromatic_atoms = list(set(aromatic_atoms))  # Remove duplicates
        
        # Build descriptive label showing which atoms are aromatic
        j_aromatic = j in aromatic_atoms
        k_aromatic = k in aromatic_atoms
        
        if j_aromatic and k_aromatic:
            # Both in ring: (Ph)C-C(Ph)
            desc = f'(Ph){symbols[j]}{j}-{symbols[k]}{k}(Ph) rotation'
        elif j_aromatic:
            # j in ring: (Ph)C-C
            desc = f'(Ph){symbols[j]}{j}-{symbols[k]}{k} rotation'
        else:
            # k in ring: C-C(Ph)
            desc = f'{symbols[j]}{j}-{symbols[k]}{k}(Ph) rotation'
        
        return {
            'type': 'aromatic_rotation',
            'description': desc,
            'axis_atoms': (j, k),
            'max_change': dihedral_change,
            'rotating_group': False,
            'aromatic_atoms': aromatic_atoms
        }
    
    def count_identical_neighbors(central: int, neighbors: List[int], exclude: int) -> Tuple[str, int]:
        """Count neighbors with same element, excluding the bond partner."""
        neighbor_symbols = [symbols[n] for n in neighbors if n != exclude]
        if not neighbor_symbols:
            return '', 0
        from collections import Counter
        counts = Counter(neighbor_symbols)
        most_common = counts.most_common(1)[0]
        return most_common
    
    # Check j end
    sym_j, count_j = count_identical_neighbors(j, neighbors_j, k)
    
    # Check k end
    sym_k, count_k = count_identical_neighbors(k, neighbors_k, j)
    
    # Identify symmetric tops (3 identical groups)
    if count_j == 3:
        rotating_atom = j
        group_symbol = sym_j
        central_symbol = symbols[j]
    elif count_k == 3:
        rotating_atom = k
        group_symbol = sym_k
        central_symbol = symbols[k]
    else:
        rotating_atom = None
        group_symbol = None
        central_symbol = None
    
    # Only classify very clear symmetric tops (3 H neighbors)
    if rotating_atom is not None and group_symbol == 'H':
        if central_symbol == 'C':
            rot_type = 'methyl'
            description = f'Methyl (CH₃) rotation at {central_symbol}{rotating_atom}'
        elif central_symbol == 'P':
            rot_type = 'phosphine'
            description = f'Phosphine (PH₃) rotation at {central_symbol}{rotating_atom}'
        elif central_symbol == 'N':
            rot_type = 'amine'
            description = f'Amine (NH₃) rotation at {central_symbol}{rotating_atom}'
        else:
            rot_type = 'symmetric_top'
            description = f'{central_symbol}H₃ rotation at {central_symbol}{rotating_atom}'
        
        return {
            'type': rot_type,
            'description': description,
            'axis_atoms': (j, k),
            'rotating_atom': rotating_atom,
            'max_change': dihedral_change,
            'rotating_group': True
        }
    
    # Generic rotation
    return {
        'type': 'single_bond',
        'description': f'Single bond {symbols[j]}-{symbols[k]} rotation',
        'axis_atoms': (j, k),
        'max_change': dihedral_change,
        'rotating_group': False
    }


def analyze_rotations(
    dihedral_changes: Dict[Tuple[int, int, int, int], Tuple[float, float]],
    frame_ts: Atoms,
    connectivity: Dict[int, set]
) -> Dict[Tuple[int, int, int, int], Dict[str, Any]]:
    """
    Analyze dihedrals and classify rotation types.
    
    Args:
        dihedral_changes: Dict of dihedrals with significant changes
        frame_ts: TS frame
        connectivity: Connectivity dictionary from build_internal_coordinates
    
    Returns:
        Dict mapping dihedrals to rotation characterizations
    """
    rotations = {}
    
    for dihedral, (max_change, initial_value) in dihedral_changes.items():
        rotation_info = classify_rotation_type(dihedral, frame_ts, max_change, connectivity)
        rotation_info['initial_value'] = initial_value
        rotations[dihedral] = rotation_info
    
    return rotations


# ============================================================================
# MAIN CHARACTERIZATION
# ============================================================================

def characterize_vib_mode(
    internal_changes: Dict[str, Any],
    frames: List[Atoms],
    ts_frame_idx: int = 0
) -> Dict[str, Any]:
    """
    Comprehensive characterization using simple, robust logic.
    
    Strategy:
    1. If bond changes detected -> bond formation/breaking
    2. If no bonds but hub pattern in dihedrals -> inversion
    3. Otherwise -> rotation(s)
    
    Args:
        internal_changes: Results from analyze_internal_displacements
        frames: Full trajectory frames
        ts_frame_idx: Index of TS frame
    
    Returns:
        Dictionary with characterization results
    """
    bond_changes = internal_changes['bond_changes']
    angle_changes = internal_changes['angle_changes']
    dihedral_changes = internal_changes['dihedral_changes']
    frame_indices = internal_changes['frame_indices']
    atom_index_map = internal_changes['atom_index_map']
    
    # Get TS and displaced frames
    frame_ts = frames[ts_frame_idx]
    frames_displaced = [frames[i] for i in frame_indices]
    symbols = frame_ts.get_chemical_symbols()
    
    # Determine primary motion type
    has_bonds = len(bond_changes) > 0
    has_angles = len(angle_changes) > 0
    has_dihedrals = len(dihedral_changes) > 0
    
    result = {
        'has_bond_changes': has_bonds,
        'has_angle_changes': has_angles,
        'has_dihedral_changes': has_dihedrals,
        'rotations': {},
        'inversion': None
    }
    
    # If bond changes, that's the primary motion
    if has_bonds:
        result['mode_type'] = 'bond_change'
        result['description'] = 'Bond formation/breaking'
        result['primary_motion'] = 'bond_change'
        logger.info(f"Characterized as bond change ({len(bond_changes)} bonds)")
        return result
    
    # No bond changes - check for inversion or rotation
    if has_dihedrals:
        # Try to detect inversion hub
        inversion_result = detect_inversion_hub(dihedral_changes, symbols)
        
        if inversion_result:
            hub_atom, element, hub_fraction = inversion_result
            
            # Identify which group is moving (using connectivity from internal_changes)
            connectivity = internal_changes.get('connectivity', {})
            moving_group = identify_moving_group(hub_atom, frame_ts, frames_displaced, connectivity)
            neighbor_idx, max_disp, group_description = moving_group
            
            result['mode_type'] = 'inversion'
            
            # Build description based on whether moving group was identified
            if group_description == "unknown" or max_disp < 0.01:
                result['description'] = f'{element} inversion'
                logger.info(f"Characterized as {element} inversion (hub atom {hub_atom})")
            else:
                result['description'] = f'{element} inversion with {group_description} motion'
                logger.info(
                    f"Characterized as {element} inversion "
                    f"(hub atom {hub_atom}, moving {group_description})"
                )
            
            result['primary_motion'] = 'inversion'
            result['inversion'] = {
                'center_atom': hub_atom,
                'center_symbol': element,
                'hub_fraction': hub_fraction,
                'moving_group': group_description if group_description != "unknown" else None,
                'moving_atom': neighbor_idx if group_description != "unknown" else None,
                'max_displacement': max_disp,
                'num_dihedrals': len(dihedral_changes)
            }
            
            return result
        
        # No inversion detected - classify as rotation(s)
        connectivity = internal_changes.get('connectivity', {})
        rotations = analyze_rotations(dihedral_changes, frame_ts, connectivity)
        result['rotations'] = rotations
        result['mode_type'] = 'rotation'
        
        # Generate description from rotation types
        rot_types = [info['type'] for info in rotations.values()]
        unique_types = list(set(rot_types))
        
        if 'methyl' in unique_types:
            result['description'] = 'Methyl rotation'
        elif 'tert_butyl' in unique_types:
            result['description'] = 'tert-Butyl rotation'
        elif 'phosphine' in unique_types:
            result['description'] = 'Phosphine rotation'
        elif len(unique_types) == 1 and unique_types[0] == 'single_bond':
            result['description'] = 'Single bond rotation'
        else:
            result['description'] = f'Multiple rotations ({", ".join(unique_types)})'
        
        result['primary_motion'] = 'rotation'
        logger.info(f"Characterized as rotation ({len(rotations)} dihedrals)")
        return result
    
    # Only angle changes (no bonds or dihedrals)
    if has_angles:
        result['mode_type'] = 'bending'
        result['description'] = 'Angle bending motion'
        result['primary_motion'] = 'bending'
        logger.info(f"Characterized as bending ({len(angle_changes)} angles)")
        return result
    
    # No significant changes detected
    result['mode_type'] = 'other'
    result['description'] = 'Other vibrational motion'
    result['primary_motion'] = 'other'
    logger.info("No significant changes detected")
    return result
