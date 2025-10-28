"""
High-level API for vibrational trajectory analysis.

This module provides the main entry points for running complete analyses,
including trajectory loading, internal coordinate analysis, and optional
graph-based analysis.
"""

import os
import logging
from typing import List, Optional, Dict, Any

from . import config
from .core import analyze_internal_displacements, read_xyz_trajectory
from .convert import (
    parse_cclib_output, 
    convert_orca, 
    get_orca_pltvib_path,
    parse_xyz_string_to_frames
)
from .graph_compare import analyze_displacement_graphs
from .utils import write_trajectory_file, save_displacement_pair, setup_logging
from .characterize import characterize_vib_mode

logger = logging.getLogger("vib_analysis")


def load_trajectory(
    input_file: str,
    mode: int = 0,
    orca_pltvib_path: Optional[str] = None,
    save_to_disk: bool = True,
    print_output: bool = False
) -> Dict[str, Any]:
    """
    Load vibrational trajectory from XYZ or QM output file.
    
    Args:
        input_file: Path to XYZ trajectory or QM output file
        mode: Vibrational mode index (ignored for XYZ files)
        orca_pltvib_path: Optional path to orca_pltvib executable
        save_to_disk: Whether to save converted trajectory to disk
        print_output: Print status messages
    
    Returns:
        Dictionary with keys:
            - 'frames': List of ASE Atoms objects
            - 'frequencies': List of frequencies (None for XYZ input)
            - 'trajectory_file': Path to trajectory file (None if not saved)
    """
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
        return {
            'frames': frames,
            'frequencies': None,
            'trajectory_file': trajectory_file
        }
    
    # QM output file - try cclib first, then ORCA
    try:
        if print_output:
            print(f"\nParsing {basename} with cclib...")
        frequencies, trajectory_string = parse_cclib_output(input_file, mode)
    except Exception as e:
        if print_output:
            print(f"cclib failed ({e}), trying orca_pltvib...")
        
        if orca_pltvib_path is None:
            orca_pltvib_path = get_orca_pltvib_path()
        
        frequencies, trajectory_string = convert_orca(
            input_file, 
            mode, 
            pltvib_path=orca_pltvib_path
        )
    
    # Convert string to frames
    frames = parse_xyz_string_to_frames(trajectory_string)
    
    # Optionally save to disk
    if save_to_disk:
        output_path = f"{root}.v{mode:03d}.xyz"
        trajectory_file = write_trajectory_file(trajectory_string, output_path)
        if print_output:
            print(f"Saved trajectory to {os.path.basename(trajectory_file)}")
    
    return {
        'frames': frames,
        'frequencies': frequencies,
        'trajectory_file': trajectory_file
    }


def run_vib_analysis(
    input_file: str,
    mode: int = 0,
    ts_frame: int = config.DEFAULT_TS_FRAME,
    # Vibrational analysis parameters
    bond_tolerance: float = config.BOND_TOLERANCE,
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    coupled_motion_filter: float = config.COUPLED_MOTION_FILTER,
    # Graph analysis parameters (includes mode characterization)
    enable_graph: bool = False,
    graph_method: str = "cheminf",
    charge: int = 0,
    multiplicity: Optional[int] = None,
    distance_tolerance: float = config.DISTANCE_TOLERANCE,
    independent_graphs: bool = False,
    ascii_scale: float = config.ASCII_SCALE,
    ascii_include_h: bool = config.ASCII_INCLUDE_H,
    ascii_neighbor_shells: int = config.ASCII_NEIGHBOR_SHELLS,
    # Output options
    save_trajectory: bool = config.SAVE_TRAJECTORY_DEFAULT,
    save_displacement: bool = config.SAVE_DISPLACEMENT_DEFAULT,
    displacement_scale: int = config.DEFAULT_DISPLACEMENT_LEVEL,
    orca_pltvib_path: Optional[str] = None,
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
        
        save_trajectory: Save converted trajectory to disk
        save_displacement: Save displaced structure pair
        displacement_scale: Displacement amplitude level (1-4)
        orca_pltvib_path: Path to orca_pltvib executable
        print_output: Print formatted analysis results to console
        show_all: Show all changes including minor angles/dihedrals
        debug: Enable debug output
        
    Returns:
        Dictionary with keys:
            - 'trajectory': Trajectory metadata (frames, frequencies, file path)
            - 'vibrational': Internal coordinate analysis results
            - 'characterization': Mode characterization (if enabled)
            - 'graph': Graph analysis results (if enabled)
            - 'displacement_files': Paths to saved displacement files (if enabled)
    """
    # Set up logging and print header if outputting to console
    if print_output or debug:
        # Print main header first
        print("=" * 80)
        print(" " * 30 + "VIB_ANALYSIS")
        print("=" * 80)
        
        # Set up logging (prints debug message if debug mode)
        setup_logging(debug=debug)
    
    # Load trajectory
    trajectory_data = load_trajectory(
        input_file,
        mode=mode,
        orca_pltvib_path=orca_pltvib_path,
        save_to_disk=save_trajectory,
        print_output=print_output
    )
    
    frames = trajectory_data['frames']
    
    if print_output:
        print(f"Loaded {len(frames)} frames from trajectory")
    
    # Analyze internal coordinates
    vib_results = analyze_internal_displacements(
        frames,
        ts_frame=ts_frame,
        bond_tolerance=bond_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        coupled_motion_filter=coupled_motion_filter,
    )
    
    # Check if anything was detected - if not, try with relaxed thresholds
    nothing_detected = (
        len(vib_results['bond_changes']) == 0 and
        len(vib_results['angle_changes']) == 0 and
        len(vib_results['dihedral_changes']) == 0
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
            bond_tolerance=bond_tolerance,
            bond_threshold=bond_threshold * 0.5,
            angle_threshold=angle_threshold * 0.5,
            dihedral_threshold=dihedral_threshold * 0.5,
            coupled_motion_filter=coupled_motion_filter * 0.5,
        )
        
        # Add metadata about relaxed thresholds
        vib_results['thresholds_relaxed'] = True
        vib_results['threshold_factor'] = 0.5
        
        if print_output:
            still_nothing = (
                len(vib_results['bond_changes']) == 0 and
                len(vib_results['angle_changes']) == 0 and
                len(vib_results['dihedral_changes']) == 0
            )
            if still_nothing:
                print("Still no significant changes detected (very small displacements).")
            else:
                print(f"Found changes with relaxed thresholds.")
    else:
        vib_results['thresholds_relaxed'] = False
    
    # Characterize mode type and run graph analysis (both enabled with --graph)
    characterization = None
    graph_results = None
    
    if enable_graph:
        if print_output:
            print("Characterizing vibrational mode...")
        
        characterization = characterize_vib_mode(
            internal_changes=vib_results,
            frames=frames,
            ts_frame_idx=ts_frame
        )
        
        logger.debug(
            f"Mode characterized as: {characterization['mode_type']} "
            f"({characterization['description']})"
        )
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
            mode_type = characterization.get('mode_type')
            if mode_type == 'rotation':
                # For rotations: highlight the rotating bond (j, k from dihedrals)
                for rot_info in characterization.get('rotations', {}).values():
                    j, k = rot_info['axis_atoms']
                    atoms_of_interest.update([j, k])
            elif mode_type == 'inversion':
                # For inversions: highlight hub atom and all its neighbors
                inv_info = characterization.get('inversion')
                if inv_info:
                    hub_atom = inv_info['center_atom']
                    atoms_of_interest.add(hub_atom)
                    # Add neighbors from connectivity (from internal_coords)
                    # This is populated by build_internal_coordinates via xyzgraph
                    connectivity = vib_results.get('connectivity', {})
                    neighbors = connectivity.get(hub_atom, set())
                    atoms_of_interest.update(neighbors)
        
        # Add ts_frame to internal_changes for graph analysis
        vib_results_with_ts = {**vib_results, 'ts_frame': ts_frame}
        
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
            print_output=print_output,
        )
    
    results_dict = {
        'trajectory': trajectory_data,
        'vibrational': vib_results,
        'characterization': characterization,
        'graph': graph_results,
        'displacement_files': displacement_files,
    }
    
    # Print formatted output if requested
    if print_output:
        from .output import print_analysis_results 
        print_analysis_results(
            results_dict,
            show_all=show_all or debug,  # Show all if requested or in debug mode
            mode=mode
        )
    
    return results_dict
