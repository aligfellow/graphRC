import os
from typing import List, Optional, Tuple, Dict, Any
from ase import Atoms
from ase.io import write


def write_trajectory_file(trajectory_string: str, output_path: str) -> str:
    """
    Write trajectory string to XYZ file.
    Returns the output path.
    """
    with open(output_path, 'w') as f:
        f.write(trajectory_string)
    return output_path


def write_displaced_structures(
    frames: List[Atoms],
    prefix: str,
    indices: Optional[List[int]] = None,
    ts_frame: int = 0,
    overwrite: bool = True
) -> List[str]:
    """
    Write displaced structures (forward/reverse) as XYZ files.
    
    Args:
        frames: List of ASE Atoms objects
        prefix: Output filename prefix
        indices: Optional list of frame indices [forward, reverse]
        ts_frame: TS frame index (used if indices is None)
        overwrite: Whether to overwrite existing files
    
    Behavior:
      - If indices provided: use first → {prefix}_F.xyz, second → {prefix}_R.xyz
      - If indices is None: use (ts_frame-1, ts_frame+1)
      - Negative indices allowed (Python style)
      - Single index → only _F written
    
    Returns:
        List of written file paths
    """
    written = []
    if not frames or len(frames) < 2:
        return written
    
    n = len(frames)
    
    def norm(i: int) -> Optional[int]:
        return i % n if -n <= i < n else None
    
    if indices is None:
        a = norm(ts_frame - 1)
        b = norm(ts_frame + 1)
        if a is None or b is None or a == b:
            a, b = 0, n - 1
        indices = [a, b]
    else:
        # Keep only first two valid normalized indices
        normed = []
        for raw in indices:
            ni = norm(raw)
            if ni is not None:
                normed.append(ni)
            if len(normed) == 2:
                break
        indices = normed
    
    if not indices:
        return written
    
    # First index (_F)
    f_idx = indices[0]
    f_path = f"{prefix}_F.xyz"
    if overwrite or not os.path.exists(f_path):
        write(f_path, frames[f_idx], format='xyz')
    written.append(f_path)
    
    # Optional second index (_R)
    if len(indices) > 1:
        r_idx = indices[1]
        if r_idx != f_idx:
            r_path = f"{prefix}_R.xyz"
            if overwrite or not os.path.exists(r_path):
                write(r_path, frames[r_idx], format='xyz')
            written.append(r_path)
    
    return written


def save_displacement_pair(
    frames: List[Atoms],
    ts_frame: int,
    output_prefix: str,
    level: int = 1,
    max_level: int = 4,
    verbose: bool = False,
) -> Optional[Tuple[str, str]]:
    """
    Save symmetric displaced structures (TS ± level) as XYZ files.
    
    Args:
        frames: List of trajectory frames
        ts_frame: Index of TS frame
        output_prefix: Prefix for output files
        level: Displacement level (1-4, corresponding to amplitudes ~0.2-0.8)
        max_level: Maximum allowed level
        verbose: Print status messages
    
    Returns:
        (forward_path, reverse_path) if successful, else None
    """
    n_frames = len(frames)
    
    if not (1 <= level <= max_level):
        if verbose:
            print(f"[save_displacement_pair] Invalid level {level} (must be 1–{max_level}).")
        return None
    
    minus_idx = ts_frame - level
    plus_idx = ts_frame + level
    
    if not (0 <= minus_idx < n_frames and 0 <= plus_idx < n_frames):
        if verbose:
            print(f"[save_displacement_pair] Level {level} out of range "
                  f"for TS {ts_frame} (total {n_frames} frames).")
        return None
    
    paths = write_displaced_structures(
        frames, 
        prefix=output_prefix, 
        indices=[minus_idx, plus_idx]
    )
    
    if verbose and len(paths) == 2:
        print(f"[save_displacement_pair] Saved displaced pair (±{level}): "
              f"{os.path.basename(paths[0])}, {os.path.basename(paths[1])}")
    
    return tuple(paths) if len(paths) == 2 else None


def print_vibrational_results(results: Dict[str, Any], show_all: bool = False):
    """Print formatted vibrational analysis results."""
    vib = results['vibrational']
    atom_map = vib.get('atom_index_map', {})
    
    def heading(title, width=80, fill='='):
        """Centered heading with padding."""
        inner = f" {title} "
        if len(inner) >= width:
            print(f"\n{inner}")
        else:
            pad = width - len(inner)
            left = pad // 2
            right = pad - left
            print("\n" + fill * left + inner + fill * right)
    
    def print_coordinate_section(title, data_dict, coord_type, unit):
        """Print a section of coordinate changes."""
        if not data_dict:
            return
        
        # Build entries
        entries = []
        for indices, (change, initial_value) in sorted(data_dict.items(), key=lambda x: -x[1][0]):
            idx_str = f"{coord_type} {indices}"
            if atom_map:
                sym_str = "[" + "-".join(atom_map[i] for i in indices) + "]"
            else:
                sym_str = ""
            entries.append((idx_str, sym_str, change, initial_value))
        
        # Compute column widths
        idx_width = max(len(e[0]) for e in entries)
        sym_width = max((len(e[1]) for e in entries), default=0)
        
        # Print header and entries
        heading(title)
        for idx_str, sym_str, change, initial_value in entries:
            print(f"{idx_str:<{idx_width}}  {sym_str:<{sym_width}}  "
                  f"Δ = {change:7.3f} {unit},  Initial = {initial_value:7.3f} {unit}")
    
    # Print main sections
    print_coordinate_section("Significant Bond Changes", vib['bond_changes'], "Bond", "Å")
    print_coordinate_section("Significant Angle Changes", vib['angle_changes'], "Angle", "°")
    print_coordinate_section("Significant Dihedral Changes", vib['dihedral_changes'], "Dihedral", "°")
    
    if vib['dihedral_changes'] and (vib['bond_changes'] or vib['angle_changes']):
        print("\nNote: Dihedrals may be artifacts of motion in the TS, "
              "not directly dependent on bond/angle changes.")
    
    # Print minor changes if requested
    if show_all:
        print_coordinate_section("Minor Angle Changes", vib['minor_angle_changes'], "Angle", "°")
        if vib['minor_angle_changes']:
            print("\nNote: These angles depend on other changes and may not be significant alone.")
        
        print_coordinate_section("Minor Dihedral Changes", vib['minor_dihedral_changes'], "Dihedral", "°")
        if vib['minor_dihedral_changes']:
            print("\nNote: These dihedrals depend on other changes and may not be significant alone.")


def print_frequency_info(frequencies: Optional[List[float]], mode: int):
    """Print information about vibrational frequencies."""
    if frequencies is None:
        return
    
    if mode < len(frequencies):
        freq = frequencies[mode]
        note = " (imaginary)" if freq < 0 else ""
        print(f"\nAnalyzed Mode {mode}: {freq:.2f} cm⁻¹{note}")
    
    # Show first 5 non-zero modes
    non_zero = [f for f in frequencies if abs(f) > 1e-5][:5]
    if non_zero:
        print("\nFirst 5 non-zero vibrational frequencies:")
        for i, freq in enumerate(non_zero):
            note = " (imaginary)" if freq < 0 else ""
            print(f"  Mode {i}: {freq:.2f} cm⁻¹{note}")
