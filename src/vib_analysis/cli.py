"""
Command-line interface for vibrational trajectory analysis.
"""
import argparse
import os
import sys
import logging

from . import config
from .api import run_vib_analysis
from .output import print_analysis_results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze vibrational trajectories for structural changes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    
    # Input/output
    parser.add_argument('input', help='Input file (XYZ trajectory or QM output)')
    
    # Mode selection
    parser.add_argument('--mode', '-m', type=int, default=0,
                       help='Vibrational mode to analyze (default: 0, ignored for XYZ)')
    parser.add_argument('--ts-frame', type=int, default=0,
                       help='Frame index to use as TS reference (default: 0)')
    
    # Vibrational analysis tolerances
    vib_group = parser.add_argument_group('vibrational analysis parameters')
    vib_group.add_argument('--bond-tolerance', type=float, default=config.BOND_TOLERANCE,
                          help=f'Bond detection tolerance factor (default: {config.BOND_TOLERANCE})')
    vib_group.add_argument('--bond-threshold', type=float, default=config.BOND_THRESHOLD,
                          help=f'Threshold for significant bond changes in Å (default: {config.BOND_THRESHOLD})')
    vib_group.add_argument('--angle-threshold', type=float, default=config.ANGLE_THRESHOLD,
                          help=f'Threshold for significant angle changes in degrees (default: {config.ANGLE_THRESHOLD})')
    vib_group.add_argument('--dihedral-threshold', type=float, default=config.DIHEDRAL_THRESHOLD,
                          help=f'Threshold for significant dihedral changes in degrees (default: {config.DIHEDRAL_THRESHOLD})')
    vib_group.add_argument('--coupled-motion-filter', type=float, default=config.COUPLED_MOTION_FILTER,
                          help=f'Coupled motion filter for filtering coupled changes in Å (default: {config.COUPLED_MOTION_FILTER}, advanced)')
    vib_group.add_argument('--coupled-proton-threshold', type=str, default=str(config.COUPLED_PROTON_THRESHOLD),
                          help=f'Reduced threshold for coupled proton transfers in Å (default: {config.COUPLED_PROTON_THRESHOLD}, use "false" to disable)')
    vib_group.add_argument('--all', '-a', action='store_true',
                          help='Report all changes including minor ones')
    
    # Graph analysis (includes mode characterization)
    graph_group = parser.add_argument_group('graph analysis parameters')
    graph_group.add_argument('--graph', '-g', action='store_true',
                            help='Enable graph-based analysis and mode characterization (rotations, inversions, aromatic systems)')
    graph_group.add_argument('--method', default='cheminf',
                            choices=['cheminf', 'xtb'],
                            help='Graph building method (default: cheminf)')
    graph_group.add_argument('--charge', type=int, default=0,
                            help='Molecular charge for graph building (default: 0)')
    graph_group.add_argument('--multiplicity', type=int,
                            help='Spin multiplicity (auto-detected if not specified)')
    graph_group.add_argument('--distance-tolerance', type=float, default=config.DISTANCE_TOLERANCE,
                            help=f'Tolerance for bond formation/breaking (default: {config.DISTANCE_TOLERANCE} Å)')
    graph_group.add_argument('--independent-graphs', '-ig', action='store_true',
                            help='Build molecular graphs from the displaced geometries rather than TS geometry with guided bonding (more rigorous for use with IRC or QRC displaced trajectories)')
    
    # ASCII visualization
    ascii_group = parser.add_argument_group('ASCII rendering options')
    ascii_group.add_argument('--ascii-scale', '-as', type=float, default=config.ASCII_SCALE,
                            help=f'Scale for ASCII molecular rendering (default: {config.ASCII_SCALE})')
    ascii_group.add_argument('--show-h', action='store_true',
                            help='Show hydrogen atoms in ASCII rendering')
    ascii_group.add_argument('--ascii-shells', '-ash', type=int, default=config.ASCII_NEIGHBOR_SHELLS,
                            help=f'Neighbor shells around transformation core (default: {config.ASCII_NEIGHBOR_SHELLS})')
    
    # Output options
    output_group = parser.add_argument_group('output options')
    output_group.add_argument('--save-displacement', '-sd', action='store_true',
                             help='Save displaced structure pair')
    output_group.add_argument('--displacement-scale', '-ds', type=int, default=config.DEFAULT_DISPLACEMENT_LEVEL,
                             help=f'Displacement level (1-{config.MAX_DISPLACEMENT_LEVEL}, ~0.2-0.8 amplitude) (default: {config.DEFAULT_DISPLACEMENT_LEVEL})')
    output_group.add_argument('--no-save', action='store_true',
                             help='Do not save trajectory to disk (keep in memory only)')
    output_group.add_argument('--orca-path', help='Path to ORCA executable directory')
    
    # Logging
    parser.add_argument('--debug', '-d', action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Parse coupled_proton_threshold (convert string to float or False)
    coupled_proton_threshold = args.coupled_proton_threshold
    if coupled_proton_threshold.lower() == 'false':
        coupled_proton_threshold = False
    else:
        try:
            coupled_proton_threshold = float(coupled_proton_threshold)
        except ValueError:
            print(f"Error: --coupled-proton-threshold must be a number or 'false', got '{coupled_proton_threshold}'", file=sys.stderr)
            sys.exit(1)
    
    # Run analysis
    try:
        results = run_vib_analysis(
            input_file=args.input,
            mode=args.mode,
            ts_frame=args.ts_frame,
            # Vibrational parameters
            bond_tolerance=args.bond_tolerance,
            bond_threshold=args.bond_threshold,
            angle_threshold=args.angle_threshold,
            dihedral_threshold=args.dihedral_threshold,
            coupled_motion_filter=args.coupled_motion_filter,
            coupled_proton_threshold=coupled_proton_threshold,
            # Graph parameters (includes mode characterization)
            enable_graph=args.graph,
            graph_method=args.method,
            charge=args.charge,
            multiplicity=args.multiplicity,
            distance_tolerance=args.distance_tolerance,
            independent_graphs=args.independent_graphs,
            ascii_scale=args.ascii_scale,
            ascii_include_h=args.show_h,
            ascii_neighbor_shells=args.ascii_shells,
            # Output options
            save_trajectory=not args.no_save,
            save_displacement=args.save_displacement,
            displacement_scale=args.displacement_scale,
            orca_pltvib_path=args.orca_path,
            print_output=True,  # CLI always prints output
            show_all=args.all,  # Show minor changes if requested
            debug=args.debug,
        )
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        if args.debug:
            raise
        sys.exit(1)
    
if __name__ == "__main__":
    main()
