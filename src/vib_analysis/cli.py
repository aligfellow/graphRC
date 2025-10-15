"""
Command-line interface for vibrational trajectory analysis.
"""
import argparse
import os
import sys
import logging

from .api import run_analysis
from .graph_compare import print_graph_analysis
from .utils import print_vibrational_results, print_frequency_info


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
    vib_group.add_argument('--bond-tolerance', type=float, default=1.4,
                          help='Bond detection tolerance factor (default: 1.4)')
    vib_group.add_argument('--angle-tolerance', type=float, default=1.1,
                          help='Angle detection tolerance factor (default: 1.1)')
    vib_group.add_argument('--dihedral-tolerance', type=float, default=1.0,
                          help='Dihedral detection tolerance factor (default: 1.0)')
    vib_group.add_argument('--bond-threshold', type=float, default=0.4,
                          help='Threshold for significant bond changes in Å (default: 0.4)')
    vib_group.add_argument('--angle-threshold', type=float, default=10.0,
                          help='Threshold for significant angle changes in degrees (default: 10.0)')
    vib_group.add_argument('--dihedral-threshold', type=float, default=20.0,
                          help='Threshold for significant dihedral changes in degrees (default: 20.0)')
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
    graph_group.add_argument('--distance-tolerance', type=float, default=0.2,
                            help='Tolerance for bond formation/breaking (default: 0.2 Å)')
    
    # ASCII visualization
    ascii_group = parser.add_argument_group('ASCII rendering options')
    ascii_group.add_argument('--ascii-scale', type=float, default=2.5,
                            help='Scale for ASCII molecular rendering (default: 2.5)')
    ascii_group.add_argument('--show-h', action='store_true',
                            help='Show hydrogen atoms in ASCII rendering')
    ascii_group.add_argument('--ascii-shells', type=int, default=1,
                            help='Neighbor shells around transformation core (default: 1)')
    
    # Output options
    output_group = parser.add_argument_group('output options')
    output_group.add_argument('--save-displacement', action='store_true',
                             help='Save displaced structure pair')
    output_group.add_argument('--level', type=int, default=1,
                             help='Displacement level (1-4, ~0.2-0.8 amplitude) (default: 1)')
    output_group.add_argument('--no-save', action='store_true',
                             help='Do not save trajectory to disk (keep in memory only)')
    output_group.add_argument('--orca-path', help='Path to ORCA executable directory')
    
    # Logging
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    parser.add_argument('-d', '--debug', action='store_true',
                       help='Enable debug output')
    
    args = parser.parse_args()
    
    # Set up logging
    log_level = logging.DEBUG if args.debug else (logging.INFO if args.verbose else logging.WARNING)
    logging.basicConfig(
        level=log_level,
        format='%(levelname)s: %(message)s'
    )
    
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
            verbose=args.verbose,
            debug=args.debug,
        )
    except Exception as e:
        print(f"Error during analysis: {e}", file=sys.stderr)
        if args.debug:
            raise
        sys.exit(1)
    
    # Print results
    print("\n" + "=" * 80)
    print(" " * 20 + "VIBRATIONAL TRAJECTORY ANALYSIS")
    print("=" * 80)
    
    # Frequency info
    frequencies = results['trajectory'].get('frequencies')
    print_frequency_info(frequencies, args.mode)
    
    # Graph analysis output (if enabled)
    if args.graph and results.get('graph'):
        print_graph_analysis(results['graph'], debug=args.debug)
    
    # Vibrational coordinate analysis
    print_vibrational_results(results, show_all=args.all)
    
    # Displacement file info
    if results.get('displacement_files'):
        f_path, r_path = results['displacement_files']
        print(f"\n✓ Saved displacement structures: "
              f"{os.path.basename(f_path)}, {os.path.basename(r_path)}")
    
    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
