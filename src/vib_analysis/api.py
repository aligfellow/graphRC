import os, sys, logging
from typing import List, Optional, Union
from .core import analyze_internal_displacements, read_xyz_trajectory
from .convert import (
    parse_cclib_output, get_orca_frequencies, convert_orca,
    get_orca_pltvib_path, parse_xyz_string_to_frames, write_displaced_structures
)
from .graph_compare import analyze_displacement_graphs

logger = logging.getLogger("vib_analysis.api")

def _load_trajectory(input_file: str,
                     mode: int,
                     no_save: bool,
                     orca_path: Optional[str],
                     print_output: bool):
    basename = os.path.basename(input_file)
    frames = None; freqs = None; trj_file_path = None; trj_data_str = None
    root, ext = os.path.splitext(input_file)
    if ext.lower() == ".xyz":
        if print_output:
            print(f"Reading trajectory from {basename}")
        frames = read_xyz_trajectory(input_file)
        trj_file_path = input_file
    else:
        # cclib first
        try:
            freqs, trj_file_path, trj_data_str = parse_cclib_output(input_file, mode, no_save=no_save, silent=True)
            if print_output:
                print(f"Parsing {basename} with cclib...")
        except Exception:
            # fallback orca_pltvib
            if print_output:
                print(f"Parsing {basename} with orca_pltvib...")
            if orca_path is None:
                pltvib_path = get_orca_pltvib_path()
            else:
                pltvib_path = os.path.join(os.path.dirname(orca_path), 'orca_pltvib')
            freqs = get_orca_frequencies(input_file)
            trj_file_path, trj_data_str = convert_orca(input_file, mode, pltvib_path=pltvib_path, no_save=no_save, silent=False)
    # normalize frames
    if trj_file_path and frames is None:
        frames = read_xyz_trajectory(trj_file_path)
    if frames is None and trj_data_str:
        frames = parse_xyz_string_to_frames(trj_data_str)
    return {
        "frames": frames,
        "freqs": freqs,
        "trj_path": trj_file_path,
        "trj_data": trj_data_str
    }

def _resolve_graph_frames(enable_graph: bool,
                          explicit_indices: Optional[List[int]],
                          vib_results: dict,
                          print_output: bool,
                          debug: bool):
    # (suppress normal prints; only emit when debugging)
    if not enable_graph:
        return None
    if explicit_indices:
        if print_output and debug:
            print(f"\n(Graph debug) Using user-specified graph frames: {explicit_indices}")
        return explicit_indices
    if 'frame_indices' in vib_results:
        fi = vib_results['frame_indices']
        if print_output and debug:
            print(f"\n(Graph debug) Using maximally-displaced frames from vibrational analysis: {fi}")
        return fi
    if print_output and debug:
        print("\n(Graph debug) No vibrational frame_indices found, using default: [0, -1]")
    return [0, -1]

def run_vib_analysis(
    input_file: str,
    mode: int = 0,
    bond_tolerance: float = 1.5,
    angle_tolerance: float = 1.1,
    dihedral_tolerance: float = 1.0,
    bond_threshold: float = 0.4,
    angle_threshold: float = 10.0,
    dihedral_threshold: float = 20.0,
    ts_frame: int = 0,
    report_all: bool = False,
    print_output: bool = False,
    orca_path: Optional[str] = None,
    save_displacement: bool = False,
    no_save: bool = False,
    frame_indices: Optional[List[int]] = None,
    graph_method: str = 'cheminf',
    charge: int = 0,
    multiplicity: Optional[int] = None,
    ascii_scale: float = 3.0,
    ascii_include_h: bool = True,
    ascii_neighbor_shells: int = 1,
    debug: bool = False,
    enable_graph: bool = False,
    save_displacement_amplitude: Optional[int] = None,
):
    loader = _load_trajectory(input_file, mode, no_save, orca_path, print_output)
    frames = loader['frames']; freqs = loader['freqs']
    if frames is None or len(frames) < 2:
        raise ValueError("Could not load trajectory or insufficient frames.")
    vib_results = analyze_internal_displacements(
        loader['trj_path'] if loader['trj_path'] else frames,
        bond_tolerance=bond_tolerance,
        angle_tolerance=angle_tolerance,
        dihedral_tolerance=dihedral_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        ts_frame=ts_frame,
    )
    # Graph analysis
    graph_results = None
    resolved_frames = _resolve_graph_frames(enable_graph, frame_indices, vib_results, print_output, debug)
    if enable_graph and resolved_frames:
        if print_output and debug:
            print(f"(Graph debug) Analyzing graph transformations for frames {resolved_frames}")
        graph_results = analyze_displacement_graphs(
            frames,
            resolved_frames,
            method=graph_method,
            charge=charge,
            multiplicity=multiplicity,
            internal_changes=vib_results,
            ascii_neighbor_shells=ascii_neighbor_shells,
            ascii_scale=ascii_scale,
            ascii_include_h=ascii_include_h,
            debug=debug
            # REMOVED: highlight_ts_bonds
        )
    if save_displacement:
        base = os.path.splitext(input_file)[0]
        if save_displacement_amplitude is not None:
            lvl = int(save_displacement_amplitude)
            n_frames = len(frames)
            if 1 <= lvl <= 4:
                minus_idx = ts_frame - lvl
                plus_idx  = ts_frame + lvl
                if 0 <= minus_idx < n_frames and 0 <= plus_idx < n_frames:
                    write_displaced_structures(
                        frames,
                        prefix=base,
                        indices=[minus_idx, plus_idx]
                    )
                    if print_output:
                        print(f"Saved displaced pair (level {lvl} ~ amp {0.2*lvl:.1f}) "
                              f"(indices {minus_idx},{plus_idx}) -> {base}_F.xyz / {base}_R.xyz")
                else:
                    if print_output:
                        print(f"Cannot save displacement level {lvl}: symmetric frames ({minus_idx},{plus_idx}) "
                              f"out of range for trajectory (0..{n_frames-1}). Adjust --ts-frame or level.")
            else:
                if print_output:
                    print(f"Ignoring invalid displacement amplitude {lvl} (must be 1–4).")
        else:
            # Fallback: nearest immediate neighbors if available
            left = ts_frame - 1
            right = ts_frame + 1
            if 0 <= left < len(frames) and 0 <= right < len(frames):
                write_displaced_structures(frames, prefix=base, indices=[left, right])
                if print_output:
                    print(f"Saved default neighboring displaced pair -> {base}_F.xyz / {base}_R.xyz")
            else:
                if print_output:
                    print("No valid neighboring frames to save displaced structures.")
    return {
        "frames": frames,
        "freqs": freqs,
        "vibrational": vib_results,
        "graph": graph_results,
        "config": {
            "mode": mode,
            "enable_graph": enable_graph,
            "graph_frames": resolved_frames,
            "save_level": save_displacement_amplitude
            # REMOVED: highlight_ts_bonds
        }
    }
