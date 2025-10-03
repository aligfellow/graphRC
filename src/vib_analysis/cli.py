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

def run_vib_analysis(
    input_file,
    mode=0,
    bond_tolerance=1.5,
    angle_tolerance=1.1,
    dihedral_tolerance=1.0,
    bond_threshold=0.4,
    angle_threshold=10.0,
    dihedral_threshold=20.0,
    ts_frame=0,
    report_all=False,
    print_output=False,
    orca_path=None,
    save_displacement=False,
    no_save=False,
    graph_analysis=False,
    debug_graph=False,
):
    # Auto-detect based on extension
    _, ext = os.path.splitext(input_file)
    basename = os.path.basename(input_file)
    trj_file_path = None
    trj_data_str = None
    frames = None
    freqs = None

    if ext.lower() == '.xyz':
        # Direct XYZ read
        if print_output:
            print(f"Reading trajectory from {basename}")
        trj_file_path = input_file
        try:
            frames = read_xyz_trajectory(input_file)
            if len(frames) < 2:
                raise ValueError("Expecting a trajectory with at least 2 frames.")
        except ValueError as e:
            print(f"Error: {e}")
            return None
    else:
        # Try cclib first (I have failed to fully supress CCLIB errors...)
        cclib_success = False
        try:
            # Get and disable all cclib loggers
            cclib_logger = logging.getLogger('cclib')
            original_handlers = cclib_logger.handlers[:]
            original_level = cclib_logger.level
            # Redirect stdout/stderr to devnull during parsing
            old_stdout = sys.stdout
            old_stderr = sys.stderr
            devnull = open(os.devnull, 'w')
            sys.stdout = devnull
            sys.stderr = devnull
            
            try:
                freqs, trj_file_path, trj_data_str = parse_cclib_output(input_file, mode, no_save=no_save, silent=True)
                cclib_success = True
            finally:
                # Restore everything
                sys.stdout = old_stdout
                sys.stderr = old_stderr
                devnull.close()
                cclib_logger.handlers = original_handlers
                cclib_logger.level = original_level
                cclib_logger.propagate = True
            
            if print_output:
                print(f"Parsing {basename} with cclib...")
                if trj_file_path and not no_save:
                    print(f"Written trajectory to: {trj_file_path}")
                elif no_save:
                    print("Trajectory kept in memory only (--no-save flag active).")
                print_first_5_nonzero_modes(freqs)
        except Exception as e_cclib:
            # Try orca_pltvib as fallback
            try:
                if print_output:
                    print(f"Parsing {basename} with orca_pltvib...")
                if orca_path is None:
                    pltvib_path = get_orca_pltvib_path()
                else:
                    pltvib_path = os.path.join(os.path.dirname(orca_path), 'orca_pltvib')
                freqs = get_orca_frequencies(input_file)
                trj_file_path, trj_data_str = convert_orca(input_file, mode, pltvib_path=pltvib_path, no_save=no_save, silent=False)
                if print_output:
                    print_first_5_nonzero_modes(freqs)
            except Exception as e_orca:
                print(f"Error: Could not parse {basename} with cclib or orca_pltvib.")
                return None

    # Prepare input for analysis: use file path if available, else parse in-memory string
    if trj_file_path:
        analysis_input = trj_file_path
        if frames is None:  # not yet loaded (parsed case)
            frames = read_xyz_trajectory(trj_file_path)
    elif trj_data_str:
        frames = parse_xyz_string_to_frames(trj_data_str)
        analysis_input = frames
    else:
        print("Error: Could not obtain trajectory data.")
        return None

    # Run analysis
    results = analyze_internal_displacements(
        analysis_input,
        bond_tolerance=bond_tolerance,
        angle_tolerance=angle_tolerance,
        dihedral_tolerance=dihedral_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        ts_frame=ts_frame,
    )

    # Graph analysis if requested (with cross-validation)
    if graph_analysis:
        graph_results = analyze_displacement_graphs(
            frames, 
            results['frame_indices'], 
            bond_tolerance=bond_tolerance,
            debug=debug_graph,
            internal_changes=results  # <--- Pass internal coords for cross-validation
        )
        results['graph_analysis'] = graph_results
        
        if print_output:
            print_graph_analysis(graph_results, results.get('atom_index_map'))

    # Save displaced structures if requested
    if save_displacement:
        base = os.path.splitext(input_file)[0]
        write_displaced_structures(frames, prefix=base, indices=[1, -1])

    # Print results
    if print_output:
        if freqs is not None:
            print(f"\nAnalysed vibrational trajectory (Mode {mode} with frequency {freqs[mode]:.2f} cm**-1):")
        else:
            print(f"\nAnalysed vibrational trajectory from {basename}:")
        print_analysis_results(results, argparse.Namespace(all=report_all))

    return results

def main():
    parser = argparse.ArgumentParser(description="Vibrational Mode Analysis Tool")
    parser.add_argument("input", help="Input file (XYZ trajectory, ORCA output, or Gaussian log)")
    parser.add_argument("--mode", type=int, default=0, help="Mode index to analyze (default: 0, zero-indexed)")
    parser.add_argument("--orca_path", type=str, help="Path to ORCA binary (optional)")
    parser.add_argument("--save-displacement", "-sd", action="store_true", help="Save slightly displaced structures (e.g. for tight optimisations from TS)")
    parser.add_argument("--no-save", action="store_true", help="Do not save trajectory file to disk (keep in memory only)")
    parser.add_argument("--graph", "-g", action="store_true", help="Perform graph analysis of structural transformations")
    parser.add_argument("--debug-graph", "-dg", action="store_true", help="Print detailed graph validation and debugging info")

    # Analysis parameters
    parser.add_argument("--bond_tolerance", type=float, default=1.4, help="Bond detection tolerance multiplier. Default: 1.4")
    parser.add_argument("--angle_tolerance", type=float, default=1.1, help="Angle detection tolerance multiplier. Default: 1.1")
    parser.add_argument("--dihedral_tolerance", type=float, default=1.0, help="Dihedral detection tolerance multiplier. Default: 1.0")
    parser.add_argument("--bond_threshold", type=float, default=0.4, help="Minimum bond change to report (Å). Default: 0.4")
    parser.add_argument("--angle_threshold", type=float, default=10.0, help="Minimum angle change to report (degrees). Default: 10")
    parser.add_argument("--dihedral_threshold", type=float, default=20.0, help="Minimum dihedral change to report (degrees). Default: 20")
    parser.add_argument("--ts_frame", type=int, default=0, help="Reference frame index. Default: 0")
    parser.add_argument("--all", action='store_true', help="Report all internal coordinate changes")

    args = parser.parse_args()

    run_vib_analysis(
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
        graph_analysis=args.graph,
        debug_graph=args.debug_graph,  # <--- new
    )

if __name__ == "__main__":
    main()
