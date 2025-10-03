import argparse
import os
from .core import analyze_internal_displacements, read_xyz_trajectory, calculate_bond_length, calculate_angle, calculate_dihedral
from .convert import parse_cclib_output, get_orca_frequencies, convert_orca, get_orca_pltvib_path

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

def print_first_5_nonzero_modes(freqs, args):
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
    mode=None,
    parse_cclib=False,
    parse_orca=False,
    bond_tolerance=1.5,
    angle_tolerance=1.1,
    dihedral_tolerance=1.0,
    bond_threshold=0.5,
    angle_threshold=10.0,
    dihedral_threshold=20.0,
    ts_frame=False,
    report_all=False,
    print_output=False,
    orca_path=None,
):

    if parse_cclib or parse_orca:
        if mode is None:
            raise ValueError("Mode index is required for Gaussian/ORCA conversion")

        if parse_cclib:
            freqs, trj_file = parse_cclib_output(input_file, mode)
        elif parse_orca:
            if orca_path is None:
                pltvib_path = get_orca_pltvib_path()
            else:
                pltvib_path = os.path.join(os.path.dirname(orca_path), 'orca_pltvib')
            freqs = get_orca_frequencies(input_file)
            trj_file = convert_orca(input_file, mode, pltvib_path=pltvib_path)
        trj_path = trj_file

        if print_output:
            print_first_5_nonzero_modes(freqs, argparse.Namespace(parse_orca=parse_orca))
    else:
        trj_file = input_file

    results = analyze_internal_displacements(
        trj_file,
        bond_tolerance=bond_tolerance,
        angle_tolerance=angle_tolerance,
        dihedral_tolerance=dihedral_tolerance,
        bond_threshold=bond_threshold,
        angle_threshold=angle_threshold,
        dihedral_threshold=dihedral_threshold,
        ts_frame=ts_frame,
    )

    if print_output and (parse_cclib or parse_orca):
        print(f"\nAnalysed vibrational trajectory (Mode {mode} with frequency {freqs[mode]} cm**-1):")
        print_analysis_results(results, argparse.Namespace(all=report_all))
    elif print_output:
        print(f"\nAnalysed vibrational trajectory from {trj_file}:")
        print_analysis_results(results, argparse.Namespace(all=report_all))

    return results

def main():
    parser = argparse.ArgumentParser(description="Vibrational Mode Analysis Tool")

    parser.add_argument("input", help="Input file (XYZ trajectory, ORCA output, or Gaussian log)")
    parser.add_argument("--parse_cclib", action="store_true", help="Process Gaussian/ORCA/other output file instead of XYZ trajectory: requires --mode (zero indexed)")
    parser.add_argument("--parse_orca", action="store_true", help="Parse ORCA output file instead of XYZ trajectory: requires --mode (zero indexed)")
    parser.add_argument("--mode", type=int, help="Mode index to analyze (for Gaussian/ORCA conversion)")
    parser.add_argument("--orca_path", type=str, help="Path to ORCA binary")

    # Analysis parameters
    parser.add_argument("--bond_tolerance", type=float, default=1.5, help="Bond detection tolerance multiplier. Default: 1.5")
    parser.add_argument("--angle_tolerance", type=float, default=1.1, help="Angle detection tolerance multiplier. Default: 1.1")
    parser.add_argument("--dihedral_tolerance", type=float, default=1.0, help="Dihedral detection tolerance multiplier. Default: 1.0")
    parser.add_argument("--bond_threshold", type=float, default=0.5, help="Minimum internal coordinate change to report. Default: 0.5")
    parser.add_argument("--angle_threshold", type=float, default=10.0, help="Minimum angle change in degrees to report. Default: 10")
    parser.add_argument("--dihedral_threshold", type=float, default=20.0, help="Minimum dihedral change in degrees to report. Default: 20")
    parser.add_argument("--ts_frame", action='store_true', default=0, help="TS frame for distances and angles in the TS. Default: 0 (first frame)")
    parser.add_argument("--all", action='store_true', default=False, help="Report all changes in angles and dihedrals.")

    args = parser.parse_args()

    run_vib_analysis(
            input_file=args.input,
            mode=args.mode,
            parse_cclib=args.parse_cclib,
            parse_orca=args.parse_orca,
            bond_tolerance=args.bond_tolerance,
            angle_tolerance=args.angle_tolerance,
            dihedral_tolerance=args.dihedral_tolerance,
            bond_threshold=args.bond_threshold,
            angle_threshold=args.angle_threshold,
            dihedral_threshold=args.dihedral_threshold,
            ts_frame=args.ts_frame,
            report_all=args.all,
            print_output=True,
            orca_path=args.orca_path if args.parse_orca else None,
        )

if __name__ == "__main__":
    main()
