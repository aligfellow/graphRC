"""
Command-line interface for vibrational trajectory analysis.
"""
import argparse
import os
import sys
import logging

from . import config
from .api import run_analysis
from .output import print_analysis_results


def main():
    parser = argparse.ArgumentParser(
        description='Analyze vibrational trajectories for structural changes',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze mode 0 of ORCA output
  %(prog)s calculation.out -m 0
  
  # Analyze XYZ trajectory with graph analysis
  %(prog)s trajectory.xyz --graph
  
  # Save displaced structures at amplitude level 2
  %(prog)s calculation.out -m 0 --save-displacement --level 2
  
  # Custom thresholds for bond/angle detection
  %(prog)s trajectory.xyz --bond-threshold 0.3 --angle-threshold 15.0
        """
    )
    
    # Input/output
    parser.add_argument('input', help='Input file (XYZ trajectory or QM output)')
    parser.add_argument('-o', '--output', help='Output file for results (not implemented)')
    
    # Mode selection
    parser.add_argument('-m', '--mode', type=int, default=0,
                       help='Vibrational mode to analyze (default: 0, ignored for XYZ)')
    parser.add_argument('--ts-frame', type=int, default=0,
                       help='Frame index to use as TS reference (default: 0)')
    
    # Vibrational analysis tolerances
    vib_group = parser.add_argument_group('vibrational analysis parameters')
    vib_group.add_argument('--bond-tolerance', type=float, default=config.BOND_TOLERANCE,
                          help=f'Bond detection tolerance factor (default: {config.BOND_TOLERANCE})')
    vib_group.add_argument('--angle-tolerance', type=float, default=config.ANGLE_TOLERANCE,
                          help=f'Angle detection tolerance factor (default: {config.ANGLE_TOLERANCE})')
    vib_group.add_argument('--dihedral-tolerance', type=float, default=config.DIHEDRAL_TOLERANCE,
                          help=f'Dihedral detection tolerance factor (default: {config.DIHEDRAL_TOLERANCE})')
    vib_group.add_argument('--bond-threshold', type=float, default=config.BOND_THRESHOLD,
                          help=f'Threshold for significant bond changes in Å (default: {config.BOND_THRESHOLD})')
    vib_group.add_argument('--angle-threshold', type=float, default=config.ANGLE_THRESHOLD,
                          help=f'Threshold for significant angle changes in degrees (default: {config.ANGLE_THRESHOLD})')
    vib_group.add_argument('--dihedral-threshold', type=float, default=config.DIHEDRAL_THRESHOLD,
                          help=f'Threshold for significant dihedral changes in degrees (default: {config.DIHEDRAL_THRESHOLD})')
    vib_group.add_argument('--bond-stability', type=float, default=config.BOND_STABILITY_THRESHOLD,
                          help=f'Bond stability threshold for filtering coupled changes in Å (default: {config.BOND_STABILITY_THRESHOLD}, advanced)')
    vib_group.add_argument('-a', '--all', action='store_true',
                          help='Report all changes including minor ones')
    
    # Graph analysis
    graph_group = parser.add_argument_group('graph analysis parameters')
    graph_group.add_argument('-g', '--graph', action='store_true',
                            help='Enable graph-based analysis')
    graph_group.add_argument('--method', default='cheminf',
                            choices=['cheminf', 'xtb'],
                            help='Graph building method (default: cheminf)')
    graph_group.add_argument('--charge', type=int, default=0,
                            help='Molecular charge for graph building (default: 0)')
    graph_group.add_argument('--multiplicity', type=int,
                            help='Spin multiplicity (auto-detected if not specified)')
    graph_group.add_argument('--distance-tolerance', type=float, default=config.DISTANCE_TOLERANCE,
                            help=f'Tolerance for bond formation/breaking (default: {config.DISTANCE_TOLERANCE} Å)')
    
    # ASCII visualization
    ascii_group = parser.add_argument_group('ASCII rendering options')
    ascii_group.add_argument('--ascii-scale', type=float, default=config.ASCII_SCALE,
                            help=f'Scale for ASCII molecular rendering (default: {config.ASCII_SCALE})')
    ascii_group.add_argument('--show-h', action='store_true',
                            help='Show hydrogen atoms in ASCII rendering')
    ascii_group.add_argument('--ascii-shells', type=int, default=config.ASCII_NEIGHBOR_SHELLS,
                            help=f'Neighbor shells around transformation core (default: {config.ASCII_NEIGHBOR_SHELLS})')
    
    # Output options
    output_group = parser.add_argument_group('output options')
    output_group.add_argument('--save-displacement', action='store_true',
                             help='Save displaced structure pair')
    output_group.add_argument('--level', type=int, default=config.DEFAULT_DISPLACEMENT_LEVEL,
                             help=f'Displacement level (1-{config.MAX_DISPLACEMENT_LEVEL}, ~0.2-0.8 amplitude) (default: {config.DEFAULT_DISPLACEMENT_LEVEL})')
    output_group.add_argument('--no-save', action='store_true',
                             help='Do not save trajectory to disk (keep in memory only)')
    output_group.add_argument('--orca-path', help='Path to ORCA executable directory')
    
    # Logging
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('-d', '--debug', action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    # Check input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' not found.", file=sys.stderr)
        sys.exit(1)
    
    # Run analysis
    try:
        results = run_analysis(
            input_file=args.input,
            mode=args.mode,
            ts_frame=args.ts_frame,
            # Vibrational parameters
            bond_tolerance=args.bond_tolerance,
            angle_tolerance=args.angle_tolerance,
            dihedral_tolerance=args.dihedral_tolerance,
            bond_threshold=args.bond_threshold,
            angle_threshold=args.angle_threshold,
            dihedral_threshold=args.dihedral_threshold,
            bond_stability_threshold=args.bond_stability,
            # Graph parameters
            enable_graph=args.graph,
            graph_method=args.method,
            charge=args.charge,
            multiplicity=args.multiplicity,
            distance_tolerance=args.distance_tolerance,
            ascii_scale=args.ascii_scale,
            ascii_include_h=args.show_h,
            ascii_neighbor_shells=args.ascii_shells,
            # Output options
            save_trajectory=not args.no_save,
            save_displacement=args.save_displacement,
            displacement_level=args.level,
            orca_pltvib_path=args.orca_path,
            print_output=True,  # CLI always prints output
            verbose=args.verbose,
            debug=args.debug,
        )
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        if args.debug:
            raise
        sys.exit(1)
    
    # Results are printed by run_analysis when print_output=True


if __name__ == "__main__":
    main()
