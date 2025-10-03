import numpy as np
from ase import Atoms
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.geometry import geometry
from ase.data import covalent_radii
from itertools import combinations
import os

def read_xyz_trajectory(file_path):
    """Reads an XYZ trajectory file and returns a list of ASE Atoms objects."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"Exiting: File {file_path} does not exist.") 
    
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
        raise ValueError(f"Exiting: Only one geometry found, make sure that this is a trj file with at least 2 frames.")
    return frames

def build_internal_coordinates(frame, bond_tolerance=1.5, angle_tolerance=1.1, dihedral_tolerance=1.1):
    """Builds internal coordinates (bonds, angles, dihedrals) for a given ASE Atoms frame."""
    cutoffs = natural_cutoffs(frame, mult=bond_tolerance)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(frame)

    bonds = []
    # Build bond list
    for i in range(len(frame)):
        indices, offsets = nl.get_neighbors(i)
        for j in indices:
            if j > i:
                bond = (int(i), int(j))
                bonds.append(bond)
    
    # Build connectivity graph from all bonds
    all_connectivity = {}
    for i, j in bonds:
        all_connectivity.setdefault(i, set()).add(j)
        all_connectivity.setdefault(j, set()).add(i)

    angles = []
    # Tighter neighbor list for angles and dihedrals
    angle_cutoffs = natural_cutoffs(frame, mult=angle_tolerance)
    angle_nl = NeighborList(angle_cutoffs, self_interaction=False, bothways=False)
    angle_nl.update(frame)

    angle_bonds = set()
    for i in range(len(frame)):
        indices, offsets = angle_nl.get_neighbors(i)
        for j in indices:
            if j > i:
                angle_bonds.add((int(i), int(j)))

    # Build connectivity for angles from stricter bonds
    angle_connectivity = {}
    for i, j in angle_bonds:
        angle_connectivity.setdefault(i, set()).add(j)
        angle_connectivity.setdefault(j, set()).add(i)

    # Build angles from angle connectivity
    angles = []
    for j in angle_connectivity:
        neighbors = list(angle_connectivity[j])
        for i, k in combinations(neighbors, 2):
            angles.append((int(i), int(j), int(k)))

    dihedrals = []
    # Build stricter bond list for dihedrals
    dihedral_cutoffs = natural_cutoffs(frame, mult=dihedral_tolerance)
    dihedral_nl = NeighborList(dihedral_cutoffs, self_interaction=False, bothways=True)
    dihedral_nl.update(frame)
    
    dihedral_bonds = []
    for i in range(len(frame)):
        indices, offsets = dihedral_nl.get_neighbors(i)
        for j in indices:
            if j > i:
                dihedral_bonds.append((int(i), int(j)))
    
    # Build connectivity for dihedrals from strictest bonds
    dihedral_connectivity = {}
    for i, j in dihedral_bonds:
        dihedral_connectivity.setdefault(i, set()).add(j)
        dihedral_connectivity.setdefault(j, set()).add(i)

    # Build dihedrals from dihedral connectivity
    dihedrals = []
    for b, c in dihedral_bonds:
        if b not in dihedral_connectivity or c not in dihedral_connectivity:
            continue
        
        a_neighbors = dihedral_connectivity[b] - {c}
        d_neighbors = dihedral_connectivity[c] - {b}

        for a in a_neighbors:
            for d in d_neighbors:
                if a != d:
                    dihedrals.append((int(a), int(b), int(c), int(d)))

    return {'bonds': bonds, 'angles': angles, 'dihedrals': dihedrals}

def calculate_bond_length(frame, i, j):
    return round(float(frame.get_distance(i, j)),3)


def calculate_angle(frame, i, j, k):
    return round(float(frame.get_angle(i, j, k, mic=True)),3)


def calculate_dihedral(frame, i, j, k, l):
    return round(float(frame.get_dihedral(i, j, k, l, mic=True)),3)


def calculate_internal_changes(frames, ts_frame, internal_coords, bond_threshold=0.5, angle_threshold=10.0, dihedral_threshold=20.0, bond_stability_threshold=0.2, angle_tolerance=1.05):
    """Tracks changes in internal coordinates across trajectory."""
    num_frames = len(frames)
    bond_changes = {}
    angle_changes = {}
    dihedral_changes = {}

    # Track max changes
    for i, j in internal_coords['bonds']:
        distances = [calculate_bond_length(frame, i, j) for frame in frames]
        max_change = round(max(distances) - min(distances),3)
        if abs(max_change) >= bond_threshold:
            # get the initial bond length
            initial_length = calculate_bond_length(ts_frame, i, j)
            bond_changes[(i, j)] = max_change, initial_length  

    changed_atoms = set()
    for bond in bond_changes:
        changed_atoms.update(bond)

    minor_angles = {}
    for i, j, k in internal_coords['angles']:
        bonds_in_angle = [tuple(sorted((i, j))), tuple(sorted((j, k)))] 
        if any(bond in bond_changes for bond in bonds_in_angle):
            continue
        if all(bond_changes.get(tuple(sorted(bond)), 0) < bond_stability_threshold for bond in bonds_in_angle):
            angles = [calculate_angle(frame, i, j, k) for frame in frames]
            max_change = round(max(angles) - min(angles), 3)
            if abs(max_change) >= angle_threshold:
                angle_atoms = set((i, j, k))
                initial_angle = calculate_angle(ts_frame, i, j, k)
                intersect_atoms = angle_atoms.intersection(changed_atoms)
                if intersect_atoms:
                    minor_angles[(i, j, k)] = max_change, initial_angle
                else:
                    angle_changes[(i, j, k)] = max_change, initial_angle

    for i, j, k, l in internal_coords['dihedrals']:
        bonds_in_dihedral = [(i, j), (j, k), (k, l)]
        if any(set(bond).issubset({i, j, k, l}) for bond in bond_changes):
            continue
        if all(bond_changes.get(tuple(sorted(bond)), 0) < bond_stability_threshold for bond in bonds_in_dihedral):
            dihedrals = [calculate_dihedral(frame, i, j, k, l) for frame in frames]
            # max_change = max(dihedrals) - min(dihedrals) # returns too large angles due to periodicity
            max_change = round(max([abs((d - dihedrals[0] + 180) % 360 - 180 ) for d in dihedrals]),3)  # Adjust for periodicity
            if max_change >= angle_threshold:
                dihedral_changes[(i, j, k, l)] = max_change

    masses = frames[0].get_masses()
    dihedral_groups = {}
    for (i, j, k, l), change in dihedral_changes.items():
        axis = tuple(sorted((j, k)))
        total_mass = masses[i] + masses[j] + masses[k] + masses[l]

        if axis not in dihedral_groups:
            dihedral_groups[axis] = []
        dihedral_groups[axis].append(((i, j, k, l), change, total_mass))    

    unique_dihedrals = {}
    dependent_dihedrals = {}

    for axis, dihedrals in dihedral_groups.items():
        dihedrals_sorted = sorted(dihedrals, key=lambda x: (x[2], x[1]), reverse=True)  # Get the one with the largest change
        dihedral, max_change, _ = dihedrals_sorted[0]
        if max_change >= dihedral_threshold:
            initial_dihedral = calculate_dihedral(ts_frame, *dihedral)
            dihedral_atoms = set(dihedral)
            intersect_atoms = dihedral_atoms.intersection(changed_atoms)
            if intersect_atoms:
                dependent_dihedrals[dihedral] = max_change, initial_dihedral
            else:
                unique_dihedrals[dihedral] = max_change, initial_dihedral

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
    xyz_file,
    bond_tolerance=1.2,
    angle_tolerance=1.1,
    dihedral_tolerance=1.0,
    bond_threshold=0.5,
    angle_threshold=10.0,
    dihedral_threshold=20.0,
    ts_frame=0,  # Default to first frame
):
    frames = read_xyz_trajectory(xyz_file)
    internal_coords = build_internal_coordinates(
        frame=frames[ts_frame],
        bond_tolerance=bond_tolerance,
        angle_tolerance=angle_tolerance,
        dihedral_tolerance=dihedral_tolerance,

    )
    selected_indices = select_most_diverse_frames(frames)
    selected_frames = [frames[i] for i in selected_indices]

    bond_changes, angle_changes, minor_angles, unique_dihedrals, dependent_dihedrals = calculate_internal_changes(
        frames=selected_frames,
        ts_frame=frames[ts_frame],
        internal_coords=internal_coords,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
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
    }

