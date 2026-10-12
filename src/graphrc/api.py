"""
High-level API for vibrational trajectory analysis.

This module provides the main entry points for running complete analyses,
including trajectory loading, internal coordinate analysis, and optional
graph-based analysis.
"""

import logging
import os
from typing import Any, Dict, Optional

import numpy as np

from . import __citation__, __version__

try:
    import xyzgraph

    _xyzg_version = xyzgraph.__version__
    _xyzg_citation = xyzgraph.__citation__
except (ImportError, AttributeError):
    _xyzg_version = "unknown"
    _xyzg_citation = None

from . import config
from .characterize import characterize_vib_mode
from .convert import (
    is_orca_output,
    parse_cclib_output,
    parse_orca_output,
    parse_xyz_string_to_frames,
    trajectory,
    validate_n_frames,
)
from .core import analyze_internal_displacements, read_xyz_trajectory
from .graph_compare import analyze_displacement_graphs
from .utils import save_displacement_pair, setup_logging, write_trajectory_file

logger = logging.getLogger("graphrc")


def collect_metadata(input_file: Optional[str], **params) -> Dict[str, Any]:
    """Collect runtime metadata for reproducibility. Only includes non-default parameters."""
    # Map parameters to their config defaults
    defaults = {
        "mode": 0,
        "ts_frame": config.DEFAULT_TS_FRAME,
        "bond_tolerance": config.BOND_TOLERANCE,
        "bond_threshold": config.BOND_THRESHOLD,
        "angle_threshold": config.ANGLE_THRESHOLD,
        "dihedral_threshold": config.DIHEDRAL_THRESHOLD,
        "enable_graph": False,
        "graph_method": "cheminf",
        "charge": 0,
        "independent_graphs": False,
    }

    # Only include non-default parameters
    non_default_params = {k: v for k, v in params.items() if v is not None and k in defaults and v != defaults[k]}

    return {
        "version": __version__,
        "citation": __citation__,
        "input_file": os.path.abspath(input_file) if input_file is not None else None,
        "xyzgraph_version": _xyzg_version,
        "xyzgraph_citation": _xyzg_citation,
        "parameters": non_default_params,
    }


def load_trajectory(
    input_file: str,
    mode: int = 0,
    n_frames: int = config.DEFAULT_N_FRAMES,
    save_to_disk: bool = True,
    print_output: bool = False,
) -> Dict[str, Any]:
    """
    Load vibrational trajectory from XYZ or QM output file.

    Args:
        input_file: Path to XYZ trajectory or QM output file
        mode: Vibrational mode index (ignored for XYZ files)
        n_frames: Frames in the generated trajectory; must be a positive multiple
            of 4. Only affects QM output — XYZ frame count is read from the file
            (default 20)
        save_to_disk: Whether to save converted trajectory to disk
        print_output: Print status messages

    Returns
    -------
        Dictionary with keys:
            - 'frames': List of ASE Atoms objects
            - 'frequencies': List of frequencies (None for XYZ input)
            - 'trajectory_file': Path to trajectory file (None if not saved)
    """
    # Validate up front so a bad value fails fast with a clear message, before any
    # file IO. n_frames only affects QM output, but is always checked for sanity.
    validate_n_frames(n_frames)

    basename = os.path.basename(input_file)
    root, ext = os.path.splitext(input_file)

    frames = None
    frequencies = None
    trajectory_file = None
    trajectory_string = None

    # Direct XYZ trajectory file
    if ext.lower() == ".xyz":
        if print_output:
            print(f"Reading trajectory from {basename}")
        frames = read_xyz_trajectory(input_file)
        trajectory_file = input_file
        return {"frames": frames, "frequencies": None, "trajectory_file": trajectory_file, "mode_displacement": None}

    # QM output file - ORCA files use the direct parser; everything else uses cclib
    mode_displacement = None
    if is_orca_output(input_file):
        try:
            if print_output:
                print(f"\nParsing {basename} with direct ORCA output parser...")
            frequencies, trajectory_string, mode_displacement = parse_orca_output(input_file, mode, n_frames)
        except Exception as e:
            if print_output:
                print(f"Direct ORCA parsing failed ({e}), falling back to cclib...")
            frequencies, trajectory_string, mode_displacement = parse_cclib_output(input_file, mode, n_frames)
    else:
        if print_output:
            print(f"\nParsing {basename} with cclib...")
        frequencies, trajectory_string, mode_displacement = parse_cclib_output(input_file, mode, n_frames)

    # Convert string to frames
    frames = parse_xyz_string_to_frames(trajectory_string)

    # Optionally save to disk
    if save_to_disk:
        output_path = f"{root}.v{mode:03d}.xyz"
        trajectory_file = write_trajectory_file(trajectory_string, output_path)
        if print_output:
            print(f"Saved trajectory to {os.path.basename(trajectory_file)}")

    return {
        "frames": frames,
        "frequencies": frequencies,
        "trajectory_file": trajectory_file,
        "mode_displacement": mode_displacement,
    }


def run_vib_analysis(
    input_file: str,
    mode: int = 0,
    ts_frame: int = config.DEFAULT_TS_FRAME,
    frame_selection: str = config.FRAME_SELECTION,
    # Vibrational analysis parameters
    relaxed: bool = config.RELAXED,
    bond_tolerance: float = config.BOND_TOLERANCE,
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    coupled_motion_filter: float = config.COUPLED_MOTION_FILTER,
    coupled_proton_threshold=config.COUPLED_PROTON_THRESHOLD,
    # Graph analysis parameters (includes mode characterization)
    enable_graph: bool = False,
    graph_method: str = "cheminf",
    charge: int = 0,
    multiplicity: Optional[int] = None,
    distance_tolerance: float = config.DISTANCE_TOLERANCE,
    independent_graphs: bool = False,
    ig_flexible: bool = config.IG_FLEXIBLE_DEFAULT,
    ascii_scale: float = config.ASCII_SCALE,
    ascii_include_h: bool = config.ASCII_INCLUDE_H,
    ascii_neighbor_shells: int = config.ASCII_NEIGHBOR_SHELLS,
    # Trajectory frame count (QM output files only; ignored for XYZ input)
    n_frames: int = config.DEFAULT_N_FRAMES,
    # Output options
    save_trajectory: bool = config.SAVE_TRAJECTORY_DEFAULT,
    save_displacement: bool = config.SAVE_DISPLACEMENT_DEFAULT,
    displacement_scale: int = config.DEFAULT_DISPLACEMENT_LEVEL,
    print_output: bool = False,
    show_all: bool = False,
    debug: bool = False,
) -> Dict[str, Any]:
    """
    Complete vibrational trajectory analysis pipeline.

    Args:
        input_file: XYZ trajectory or QM output file
        mode: Vibrational mode to analyze
        ts_frame: Frame index to use as TS reference
        frame_selection: Frame selection method ('rmsd' or 'bookend')

        bond_tolerance: Multiplier for bond detection cutoffs
        bond_threshold: Threshold for significant bond changes (Å)
        angle_threshold: Threshold for significant angle changes (degrees)
        dihedral_threshold: Threshold for significant dihedral changes (degrees)
        coupled_motion_filter: Threshold for filtering coupled angle/dihedral changes (Å)

        enable_characterization: Enable mode characterization (rotations, inversions)

        enable_graph: Enable graph-based analysis
        graph_method: Method for graph building ('cheminf' or 'xtb')
        charge: Molecular charge for graph building
        multiplicity: Spin multiplicity (auto-detected if None)
        distance_tolerance: Tolerance for bond formation/breaking detection (Å)
        independent_graphs: Build graphs from actual displaced geometries rather than TS geometry
        ascii_scale: Scale factor for ASCII molecular rendering
        ascii_include_h: Include hydrogens in ASCII rendering
        ascii_neighbor_shells: Neighbor shells around transformation core

        n_frames: Frames in the generated trajectory; must be a positive multiple
            of 4. Only affects QM output — XYZ frame count is read from the file
            (default 20)
        save_trajectory: Save converted trajectory to disk
        save_displacement: Save displaced structure pair
        displacement_scale: Displacement amplitude level (1-4)
        print_output: Print formatted analysis results to console
        show_all: Show all changes including minor angles/dihedrals
        debug: Enable debug output

    Returns
    -------
        Dictionary with keys:
            - 'trajectory': Trajectory metadata (frames, frequencies, file path)
            - 'vibrational': Internal coordinate analysis results
            - 'characterization': Mode characterization (if enabled)
            - 'graph': Graph analysis results (if enabled)
            - 'displacement_files': Paths to saved displacement files (if enabled)
    """
    # Fail fast on an invalid frame count before printing anything or doing work.
    validate_n_frames(n_frames)

    if print_output or debug:
        print("=" * 80)
        print(" " * 32 + "GRAPHRC")
        print(" " * 12 + "Internal Coordinate Analysis of Vibrational Modes")
        print(" " * 26 + "A. S. Goodfellow, 2025")
        print("=" * 80)

    # Print metadata header before loading so nothing bleeds into it
    metadata = collect_metadata(
        input_file=input_file,
        mode=mode,
        ts_frame=ts_frame,
        bond_tolerance=bond_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        enable_graph=enable_graph,
        graph_method=graph_method if enable_graph else None,
        charge=charge if enable_graph else None,
        independent_graphs=independent_graphs if enable_graph else None,
    )

    if print_output:
        from .output import print_metadata_header

        print_metadata_header(metadata, {})

    if debug:
        print("Initiating debugging:")

    if print_output or debug:
        setup_logging(debug=debug)

    trajectory_data = load_trajectory(
        input_file,
        mode=mode,
        n_frames=n_frames,
        save_to_disk=save_trajectory,
        print_output=False,
    )

    frames = trajectory_data["frames"]
    mode_displacement = trajectory_data.get("mode_displacement")

    # Now print loading info after metadata will be displayed
    if print_output:
        print(f"Reading trajectory from {os.path.basename(input_file)}")
        print(f"Loaded {len(frames)} frames from trajectory")
        print(f"Using TS frame: {ts_frame}")

    analysis_results = _analyse_frames(
        frames,
        ts_frame=ts_frame,
        frame_selection=frame_selection,
        relaxed=relaxed,
        bond_tolerance=bond_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        coupled_motion_filter=coupled_motion_filter,
        coupled_proton_threshold=coupled_proton_threshold,
        enable_graph=enable_graph,
        graph_method=graph_method,
        charge=charge,
        multiplicity=multiplicity,
        distance_tolerance=distance_tolerance,
        independent_graphs=independent_graphs,
        ig_flexible=ig_flexible,
        ascii_scale=ascii_scale,
        ascii_include_h=ascii_include_h,
        ascii_neighbor_shells=ascii_neighbor_shells,
        print_output=print_output,
        debug=debug,
    )

    # Save displacement structures (optional)
    displacement_files = None
    if save_displacement:
        output_prefix = os.path.splitext(os.path.basename(input_file))[0]
        displacement_files = save_displacement_pair(
            frames=frames,
            ts_frame=ts_frame,
            output_prefix=output_prefix,
            scale=displacement_scale,
            max_level=config.MAX_DISPLACEMENT_LEVEL,
            mode_displacement=mode_displacement,
            print_output=print_output,
        )

    # Build results dictionary (metadata already collected earlier)
    results_dict = {
        "metadata": metadata,
        "trajectory": trajectory_data,
        **analysis_results,
        "displacement_files": displacement_files,
    }

    # Print formatted output if requested
    if print_output:
        from .output import print_analysis_results

        print_analysis_results(
            results_dict,
            show_all=show_all or debug,  # Show all if requested or in debug mode
            mode=mode,
        )

    return results_dict


def analyse(
    symbols,
    positions,
    mode,
    *,
    charge: int = 0,
    multiplicity: Optional[int] = None,
    frequency: Optional[float] = None,
    enable_graph: bool = False,
    n_frames: int = config.DEFAULT_N_FRAMES,
    **thresholds,
) -> Dict[str, Any]:
    """Analyze arrays using the same parameters and result structure as ``run_vib_analysis``.

    Positions and mode are (n_atoms, 3) arrays in Å, used without normalization.
    Frames span ``positions ± mode``; frequency only labels results in cm⁻¹.
    No trajectories or displacement pairs are saved; xTB uses scratch files.
    """
    frames = trajectory(symbols, positions, mode, n_frames)
    analysis_results = _analyse_frames(
        frames,
        charge=charge,
        multiplicity=multiplicity,
        enable_graph=enable_graph,
        print_output=False,
        debug=False,
        **thresholds,
    )
    metadata_params = {
        k: v for k, v in thresholds.items() if enable_graph or k not in ("graph_method", "independent_graphs")
    }
    return {
        "metadata": collect_metadata(
            None, enable_graph=enable_graph, charge=charge if enable_graph else None, **metadata_params
        ),
        "trajectory": {
            "frames": frames,
            "frequencies": [frequency] if frequency is not None else None,
            "trajectory_file": None,
            "mode_displacement": np.asarray(mode, dtype=float).copy(),
        },
        **analysis_results,
        "displacement_files": None,
    }


def _analyse_frames(
    frames,
    *,
    ts_frame: int = config.DEFAULT_TS_FRAME,
    frame_selection: str = config.FRAME_SELECTION,
    relaxed: bool = config.RELAXED,
    bond_tolerance: float = config.BOND_TOLERANCE,
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    coupled_motion_filter: float = config.COUPLED_MOTION_FILTER,
    coupled_proton_threshold=config.COUPLED_PROTON_THRESHOLD,
    enable_graph: bool = False,
    graph_method: str = "cheminf",
    charge: int = 0,
    multiplicity: Optional[int] = None,
    distance_tolerance: float = config.DISTANCE_TOLERANCE,
    independent_graphs: bool = False,
    ig_flexible: bool = config.IG_FLEXIBLE_DEFAULT,
    ascii_scale: float = config.ASCII_SCALE,
    ascii_include_h: bool = config.ASCII_INCLUDE_H,
    ascii_neighbor_shells: int = config.ASCII_NEIGHBOR_SHELLS,
    print_output: bool = False,
    debug: bool = False,
) -> Dict[str, Any]:
    """Run the shared analysis pipeline on frames, without file loading or saving."""
    # Analyze internal coordinates
    vib_results = analyze_internal_displacements(
        frames,
        ts_frame=ts_frame,
        frame_selection=frame_selection,
        relaxed=relaxed,
        bond_tolerance=bond_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        coupled_motion_filter=coupled_motion_filter,
        coupled_proton_threshold=coupled_proton_threshold,
        independent_graphs=independent_graphs,
        ig_flexible=ig_flexible,
    )

    if print_output:
        selected_frames = vib_results.get("frame_indices", [])
        print(f"Selected diverse frames for analysis: {selected_frames}")

    # Check if anything was detected - if not, try with relaxed thresholds
    nothing_detected = (
        len(vib_results["bond_changes"]) == 0
        and len(vib_results["angle_changes"]) == 0
        and len(vib_results["dihedral_changes"]) == 0
    )

    if nothing_detected:
        if print_output:
            print("\nNo significant changes detected with standard thresholds.")
            print("Relaxing criteria (50% thresholds)...")

        logger.info("No changes detected, retrying with 50% thresholds")

        # Retry with 50% thresholds
        vib_results = analyze_internal_displacements(
            frames,
            ts_frame=ts_frame,
            frame_selection=frame_selection,
            bond_tolerance=bond_tolerance,
            bond_threshold=bond_threshold * 0.5,
            angle_threshold=angle_threshold * 0.5,
            dihedral_threshold=dihedral_threshold * 0.5,
            coupled_motion_filter=coupled_motion_filter * 0.5,
            coupled_proton_threshold=coupled_proton_threshold,
        )

        # Add metadata about relaxed thresholds
        vib_results["thresholds_relaxed"] = True
        vib_results["threshold_factor"] = 0.5

        if print_output:
            still_nothing = (
                len(vib_results["bond_changes"]) == 0
                and len(vib_results["angle_changes"]) == 0
                and len(vib_results["dihedral_changes"]) == 0
            )
            if still_nothing:
                print("Still no significant changes detected (very small displacements).")
            else:
                print("Found changes with relaxed thresholds.")
    else:
        vib_results["thresholds_relaxed"] = False

    # Characterize mode type and run graph analysis (both enabled with --graph)
    characterization = None
    graph_results = None

    if enable_graph:
        if print_output:
            print("Characterizing vibrational mode...")

        characterization = characterize_vib_mode(internal_changes=vib_results, frames=frames, ts_frame_idx=ts_frame)

        logger.debug(f"Mode characterized as: {characterization['mode_type']} ({characterization['description']})")
        if print_output:
            print("Running graph-based analysis...")
            # Display graph building mode
            if independent_graphs:
                print("Graph building mode: Independent (from actual displaced geometries)")
            else:
                print("Graph building mode: TS-centric (with guided bonding)")

        # Extract atoms of interest from characterization for ASCII highlighting
        atoms_of_interest = set()
        if characterization:
            mode_type = characterization.get("mode_type")
            if mode_type == "rotation":
                # For rotations: highlight the rotating bond (j, k from dihedrals)
                for rot_info in characterization.get("rotations", {}).values():
                    j, k = rot_info["axis_atoms"]
                    atoms_of_interest.update([j, k])
            elif mode_type == "inversion":
                # For inversions: highlight hub atom and all its neighbors
                inv_info = characterization.get("inversion")
                if inv_info:
                    hub_atom = inv_info["center_atom"]
                    atoms_of_interest.add(hub_atom)
                    # Add neighbors from connectivity (from internal_coords)
                    # This is populated by build_internal_coordinates via xyzgraph
                    connectivity = vib_results.get("connectivity", {})
                    neighbors = connectivity.get(hub_atom, set())
                    atoms_of_interest.update(neighbors)

        # Add ts_frame to internal_changes for graph analysis
        vib_results_with_ts = {**vib_results, "ts_frame": ts_frame}

        graph_results = analyze_displacement_graphs(
            frames=frames,
            internal_changes=vib_results_with_ts,
            atoms_of_interest=list(atoms_of_interest) if atoms_of_interest else None,
            method=graph_method,
            charge=charge,
            multiplicity=multiplicity,
            distance_tolerance=distance_tolerance,
            independent_graphs=independent_graphs,
            ascii_scale=ascii_scale,
            ascii_include_h=ascii_include_h,
            ascii_neighbor_shells=ascii_neighbor_shells,
            debug=debug,
        )

    return {
        "vibrational": vib_results,
        "characterization": characterization,
        "graph": graph_results,
    }
