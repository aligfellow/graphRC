"""
High-level API for vibrational trajectory analysis.

This module provides the main entry points for running complete analyses,
including trajectory loading, internal coordinate analysis, and optional
graph-based analysis.
"""

import os
import logging
from typing import List, Optional, Dict, Any
from ase import Atoms

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
            print(f"Parsing {basename} with cclib...")
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
    angle_tolerance: float = config.ANGLE_TOLERANCE,
    dihedral_tolerance: float = config.DIHEDRAL_TOLERANCE,
    bond_threshold: float = config.BOND_THRESHOLD,
    angle_threshold: float = config.ANGLE_THRESHOLD,
    dihedral_threshold: float = config.DIHEDRAL_THRESHOLD,
    bond_stability_threshold: float = config.BOND_STABILITY_THRESHOLD,
    # Graph analysis parameters
    enable_graph: bool = False,
    graph_method: str = "cheminf",
    charge: int = 0,
    multiplicity: Optional[int] = None,
    distance_tolerance: float = config.DISTANCE_TOLERANCE,
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
        angle_tolerance: Multiplier for angle detection cutoffs
        dihedral_tolerance: Multiplier for dihedral detection cutoffs
        bond_threshold: Threshold for significant bond changes (Å)
        angle_threshold: Threshold for significant angle changes (degrees)
        dihedral_threshold: Threshold for significant dihedral changes (degrees)
        bond_stability_threshold: Threshold for filtering coupled angle/dihedral changes (Å)
        
        enable_graph: Enable graph-based analysis
        graph_method: Method for graph building ('cheminf' or 'xtb')
        charge: Molecular charge for graph building
        multiplicity: Spin multiplicity (auto-detected if None)
        distance_tolerance: Tolerance for bond formation/breaking detection (Å)
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
        angle_tolerance=angle_tolerance,
        dihedral_tolerance=dihedral_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        bond_stability_threshold=bond_stability_threshold,
    )
    
    # Graph analysis (optional)
    graph_results = None
    if enable_graph:
        if print_output:
            print("Running graph-based analysis...")
        
        # Add ts_frame to internal_changes for graph analysis
        vib_results_with_ts = {**vib_results, 'ts_frame': ts_frame}
        
        graph_results = analyze_displacement_graphs(
            frames=frames,
            internal_changes=vib_results_with_ts,
            method=graph_method,
            charge=charge,
            multiplicity=multiplicity,
            distance_tolerance=distance_tolerance,
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
