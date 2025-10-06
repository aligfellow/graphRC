import argparse
import os
import sys
import tempfile
import logging
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from .core import analyze_internal_displacements, read_xyz_trajectory, calculate_bond_length, calculate_angle, calculate_dihedral
from .convert import parse_cclib_output, get_orca_frequencies, convert_orca, get_orca_pltvib_path, write_displaced_structures, parse_xyz_string_to_frames
from .graph_compare import analyze_displacement_graphs, print_graph_analysis
from .api import run_vib_analysis  # NEW import (moved logic)

def print_analysis_results(results, args):
    atom_map = results.get('atom_index_map')

    # --- added helper for consistent-width headings ---
    def heading(title, total=80, fill='='):
        inner = f" {title} "
        if len(inner) >= total:
            print(f"\n{inner}")
        else:
            pad = total - len(inner)
            left = pad // 2
            right = pad - left
            print("\n" + fill * left + inner + fill * right)

    def section(title, data_dict, kind, unit):
        if not data_dict:
            return
        # Build entries first
        entries = []
        for indices, (change, initial_value) in sorted(data_dict.items(), key=lambda x: -x[1][0]):
            idx_str = f"{kind} {indices}"
            if atom_map:
                sym_str = "[" + "-".join(atom_map[i] for i in indices) + "]"
            else:
                sym_str = ""
            entries.append((idx_str, sym_str, change, initial_value))
        # Compute widths
        idx_w = max(len(e[0]) for e in entries)
        sym_w = max((len(e[1]) for e in entries), default=0)
        # Header
        heading(title)
        # Lines
        for idx_str, sym_str, change, initial_value in entries:
            print(f"{idx_str:<{idx_w}}  {sym_str:<{sym_w}}  Δ = {change:7.3f} {unit},  Initial = {initial_value:7.3f} {unit}")

    section("Significant Bond Changes", results['bond_changes'], "Bond", "Å")
    section("Significant Angle Changes", results['angle_changes'], "Angle", "°")
    section("Significant Dihedral Changes", results['dihedral_changes'], "Dihedral", "°")

    if results['dihedral_changes'] and (results['bond_changes'] or results['angle_changes']):
        print("\nNote: These dihedrals are not directly dependent on other changes however they may be artefacts of motion in the TS.")

    if args.all:
        section("Minor Angle Changes", results['minor_angle_changes'], "Angle", "°")
        if results['minor_angle_changes']:
            print("\nNote: These angles are dependent on other changes and may not be significant on their own.")
        section("Minor Dihedral Changes", results['minor_dihedral_changes'], "Dihedral", "°")
        if results['minor_dihedral_changes']:
            print("\nNote: These dihedrals are dependent on other changes and may not be significant on their own.")

def print_first_5_nonzero_modes(freqs):
    """Print the first 5 non-zero vibrational modes with proper handling"""
    # Filter out zero frequencies and get first 5
    non_zero = [f for f in freqs if abs(f) > 1e-5][:5]
    
    print("\nFirst 5 non-zero vibrational frequencies:")
    for i, freq in enumerate(non_zero):
        # Add note for imaginary frequencies
        note = " (imaginary)" if freq < 0 else ""
        print(f"  Mode {i}: {freq:.2f} cm**-1 {note}")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze vibrational trajectories for bond changes'
    )
    
    # Input/output
    parser.add_argument('input', help='Input file (XYZ trajectory or QM output)')
    parser.add_argument('-o', '--output', help='Output file for results')
    
    # Vibrational analysis parameters (for internal coordinates)
    parser.add_argument('-m', '--mode', type=int, default=0,
                       help='Vibrational mode to analyze (default: 0)')
    parser.add_argument('--bond-tolerance', type=float, default=1.4,
                       help='Bond detection tolerance factor (default: 1.4)')
    parser.add_argument('--angle-tolerance', type=float, default=1.1,
                       help='Angle detection tolerance factor (default: 1.1)')
    parser.add_argument('--dihedral-tolerance', type=float, default=1.0,
                       help='Dihedral detection tolerance factor (default: 1.0)')
    parser.add_argument('--bond-threshold', type=float, default=0.4,
                       help='Threshold for significant bond changes in Å (default: 0.4)')
    parser.add_argument('--angle-threshold', type=float, default=10.0,
                       help='Threshold for significant angle changes in degrees (default: 10.0)')
    parser.add_argument('--dihedral-threshold', type=float, default=20.0,
                       help='Threshold for significant dihedral changes in degrees (default: 20.0)')
    parser.add_argument('--ts-frame', type=int, default=0,
                       help='Frame index to use as TS reference (default: 0)')
    parser.add_argument('-a', '--all', action='store_true',
                       help='Report all changes including minor ones')
    
    # Graph analysis parameters (separate from vibrational)
    parser.add_argument('-g', '--graph', dest='enable_graph', action='store_true',
                       help='Enable graph comparison (uses maximally-displaced frames from vib analysis)')
    parser.add_argument('-f', '--frames', type=int, nargs='+',
                       help='Specific frame indices for graph comparison (overrides auto-detection)')
    parser.add_argument('--method', default='cheminf',
                       choices=['cheminf', 'xtb'],
                       help='Graph building method (default: cheminf)')
    parser.add_argument('--charge', type=int, default=0,
                       help='Molecular charge for graph building (default: 0)')
    parser.add_argument('--multiplicity', type=int,
                       help='Spin multiplicity (auto-detect if not specified)')
    
    # ASCII rendering options
    parser.add_argument('--ascii-scale', type=float, default=2.5,
                       help='Scale for ASCII molecular rendering (default: 3.0)')
    parser.add_argument('--show-h', action='store_true',
                       help='Show hydrogen atoms in ASCII rendering')
    parser.add_argument('--ascii-shells', type=int, default=1,
                       help='Neighbor shells around transformation core (default: 1)')
    
    # Other options
    parser.add_argument('--save-displacement', action='store_true',
                       help='Save displaced structures (default nearest forward/reverse or specified level)')
    parser.add_argument('--displacement-amplitude', type=int,
                       help='Symmetric displacement level (1–4 => amplitudes ~0.2–0.8) about --ts-frame; saves two frames.')
    parser.add_argument('--no-save', action='store_true',
                       help='Keep trajectory in memory only (don\'t save to disk)')
    parser.add_argument('--orca-path', help='Path to ORCA executable (for orca_pltvib)')
    parser.add_argument('--debug', '-d', action='store_true',
                       help='Enable debug output for graph analysis')
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    
    args = parser.parse_args()
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(level=log_level, format='%(levelname)s: %(message)s')
    frame_indices = args.frames if args.frames else None

    results = run_vib_analysis(
        input_file=args.input,
        mode=args.mode,
        bond_tolerance=args.bond_tolerance,
        angle_tolerance=args.angle_tolerance,
        dihedral_tolerance=args.dihedral_tolerance,
        bond_threshold=args.bond_threshold,
        angle_threshold=args.angle_threshold,
        dihedral_threshold=args.dihedral_threshold,
        ts_frame=args.ts_frame,
        report_all=args.all,
        print_output=True,
        orca_path=args.orca_path,
        save_displacement=args.save_displacement,
        no_save=args.no_save,
        frame_indices=frame_indices,
        graph_method=args.method,
        charge=args.charge,
        multiplicity=args.multiplicity,
        ascii_scale=args.ascii_scale,
        ascii_include_h=args.show_h,
        ascii_neighbor_shells=args.ascii_shells,
        debug=args.debug,
        enable_graph=args.enable_graph,
        save_displacement_amplitude=args.displacement_amplitude, 
    )

    vib = results['vibrational']
    freqs = results['freqs']

    # Print graph analysis first (if enabled)
    if args.enable_graph and results.get('graph'):
        from .graph_compare import print_graph_analysis
        print_graph_analysis(results['graph'], debug=args.debug)

    # Now print vibrational analysis so significant changes appear last
    if freqs is not None:
        print(f"\nAnalysed vibrational trajectory (Mode {args.mode} with frequency {freqs[args.mode]:.2f} cm**-1):")
    else:
        print(f"\nAnalysed vibrational trajectory from {os.path.basename(args.input)}:")
    print_analysis_results(vib, argparse.Namespace(all=args.all))

if __name__ == "__main__":
    main()
