"""
Output formatting and printing functions for vibrational analysis.

This module handles all formatted output for CLI and interactive use,
including vibrational coordinate changes, graph analysis results, and
frequency information.
"""

import os
from typing import Dict, Any, Optional, List
import logging

logger = logging.getLogger("vib_analysis")


# =====================================================================================
# === MASTER OUTPUT FUNCTION ===
# =====================================================================================

def print_analysis_results(
    results: Dict[str, Any],
    show_all: bool = False,
    mode: int = 0
) -> None:
    """
    Print complete analysis results with proper formatting.
    
    Centralized output function used by both CLI and API.
    Provides consistent formatting for all analysis results.
    Note: VIB_ANALYSIS header is printed by run_analysis() before this function.
    
    Args:
        results: Results dictionary from run_analysis()
        show_all: If True, include minor angle and dihedral changes
        mode: Vibrational mode number (for frequency display)
    """
    # Frequency info
    frequencies = results['trajectory'].get('frequencies')
    print_frequency_info(frequencies, mode)
    
    # Graph analysis (if present)
    if results.get('graph'):
        print_graph_analysis(results['graph'])
    
    # Vibrational trajectory header
    print("\n" + "=" * 80)
    print(" " * 20 + "VIBRATIONAL TRAJECTORY ANALYSIS")
    print("=" * 80)
    
    # Vibrational results
    print_vibrational_results(results, show_all=show_all)
    
    # Displacement file info (if present)
    if results.get('displacement_files'):
        f_path, r_path = results['displacement_files']
        print(f"\n✓ Saved displacement structures: "
              f"{os.path.basename(f_path)}, {os.path.basename(r_path)}")
    
    print("\n" + "=" * 80)


# =====================================================================================
# === VIBRATIONAL ANALYSIS OUTPUT ===
# =====================================================================================

def print_vibrational_results(results: Dict[str, Any], show_all: bool = False) -> None:
    """
    Print formatted vibrational analysis results.
    
    Args:
        results: Results dictionary from run_analysis containing 'vibrational' key
        show_all: If True, include minor angle and dihedral changes
    """
    vib = results['vibrational']
    atom_map = vib.get('atom_index_map', {})
    
    def heading(title: str, width: int = 80, fill: str = '=') -> None:
        """Print centered heading with padding."""
        inner = f" {title} "
        if len(inner) >= width:
            print(f"\n{inner}")
        else:
            pad = width - len(inner)
            left = pad // 2
            right = pad - left
            print("\n" + fill * left + inner + fill * right)
    
    def print_coordinate_section(
        title: str, 
        data_dict: Dict, 
        coord_type: str, 
        unit: str
    ) -> None:
        """
        Print a section of coordinate changes.
        
        Args:
            title: Section title
            data_dict: Dictionary of coordinate changes
            coord_type: Type name (e.g., "Bond", "Angle", "Dihedral")
            unit: Unit string (e.g., "Å", "°")
        """
        if not data_dict:
            return
        
        # Build entries
        entries = []
        for indices, (change, initial_value) in sorted(
            data_dict.items(), 
            key=lambda x: -x[1][0]
        ):
            idx_str = f"{coord_type} {indices}"
            if atom_map:
                sym_str = "[" + "-".join(atom_map[i] for i in indices) + "]"
            else:
                sym_str = ""
            entries.append((idx_str, sym_str, change, initial_value))
        
        # Compute column widths
        idx_width = max(len(e[0]) for e in entries)
        sym_width = max((len(e[1]) for e in entries), default=0)
        
        # Print header and entries
        heading(title)
        for idx_str, sym_str, change, initial_value in entries:
            print(
                f"{idx_str:<{idx_width}}  {sym_str:<{sym_width}}  "
                f"Δ = {change:7.3f} {unit},  Initial = {initial_value:7.3f} {unit}"
            )
    
    # Print main sections
    print_coordinate_section(
        "Significant Bond Changes", 
        vib['bond_changes'], 
        "Bond", 
        "Å"
    )
    print_coordinate_section(
        "Significant Angle Changes", 
        vib['angle_changes'], 
        "Angle", 
        "°"
    )
    print_coordinate_section(
        "Significant Dihedral Changes", 
        vib['dihedral_changes'], 
        "Dihedral", 
        "°"
    )
    
    if vib['dihedral_changes'] and (vib['bond_changes'] or vib['angle_changes']):
        print(
            "\nNote: Dihedrals may be artifacts of motion in the TS, "
            "not directly dependent on bond/angle changes."
        )
    
    # Print minor changes if requested
    if show_all:
        print_coordinate_section(
            "Minor Angle Changes", 
            vib['minor_angle_changes'], 
            "Angle", 
            "°"
        )
        if vib['minor_angle_changes']:
            print(
                "\nNote: These angles depend on other changes and may not be "
                "significant alone."
            )
        
        print_coordinate_section(
            "Minor Dihedral Changes", 
            vib['minor_dihedral_changes'], 
            "Dihedral", 
            "°"
        )
        if vib['minor_dihedral_changes']:
            print(
                "\nNote: These dihedrals depend on other changes and may not be "
                "significant alone."
            )


def print_frequency_info(frequencies: Optional[List[float]], mode: int) -> None:
    """
    Print information about vibrational frequencies.
    
    Args:
        frequencies: List of frequencies from QM calculation (None for XYZ input)
        mode: Mode index that was analyzed
    """
    if frequencies is None:
        return
    
    if mode < len(frequencies):
        freq = frequencies[mode]
        note = " (imaginary)" if freq < 0 else ""
        print(f"\nAnalyzed Mode {mode}: {freq:.2f} cm⁻¹{note}")
    
    # Show first 5 non-zero modes
    non_zero = [f for f in frequencies if abs(f) > 1e-5][:5]
    if non_zero:
        print("\nFirst 5 non-zero vibrational frequencies:")
        for i, freq in enumerate(non_zero):
            note = " (imaginary)" if freq < 0 else ""
            print(f"  Mode {i}: {freq:.2f} cm⁻¹{note}")


# =====================================================================================
# === GRAPH ANALYSIS OUTPUT ===
# =====================================================================================

def interpret_bond_order(order: float) -> str:
    """
    Convert numeric bond order to readable string.
    
    Args:
        order: Numeric bond order value
        
    Returns:
        Human-readable bond order description
    """
    if abs(order - 1.0) < 0.1:
        return "single"
    elif abs(order - 1.5) < 0.1:
        return "aromatic"
    elif abs(order - 2.0) < 0.1:
        return "double"
    elif abs(order - 3.0) < 0.1:
        return "triple"
    else:
        return f"{order:.1f}"


def print_graph_analysis(
    results: Dict[str, Any],
    atom_index_map: Optional[Dict[int, str]] = None,
    debug: bool = False
) -> None:
    """
    Print formatted graph comparison summary with optional ASCII output.
    
    Args:
        results: Results dictionary from analyze_displacement_graphs
        atom_index_map: Optional mapping of atom indices to symbols
        debug: If True, include ASCII visualization
    """
    comp = results["comparison"]
    g_ts = results["ts_graph"]
    g1 = results["frame1_graph"]
    g2 = results["frame2_graph"]

    print("\n" + "=" * 80)
    print(" " * 25 + "VIBRATIONAL GRAPH ANALYSIS SUMMARY")
    print("=" * 80)

    print("\nInterpret with care, bond orders and charges are empirical and LOW confidence.")

    # Bond formation / breaking
    formed = comp.get("bonds_formed", [])
    broken = comp.get("bonds_broken", [])
    bond_order_changes = comp.get("bond_order_changes", {})
    charge_shift = comp.get("charge_redistribution", {})

    if formed:
        print(f"\nBonds Formed ({len(formed)}):")
        for (i, j) in formed:
            s1 = g2.nodes[i].get("symbol", "?")
            s2 = g2.nodes[j].get("symbol", "?")
            order2 = g2[i][j].get('bond_order', 1)
            order_str = interpret_bond_order(order2)
            print(
                f"  Bond ({i}, {j}) [{s1}-{s2}]: "
                f"formed as {order_str} (order={order2:.1f})"
            )

    if broken:
        print(f"\nBonds Broken ({len(broken)}):")
        for (i, j) in broken:
            s1 = g1.nodes[i].get("symbol", "?")
            s2 = g1.nodes[j].get("symbol", "?")
            order1 = g1[i][j].get('bond_order', 1)
            order_str = interpret_bond_order(order1)
            print(
                f"  Bond ({i}, {j}) [{s1}-{s2}]: "
                f"broken from {order_str} (order={order1:.1f})"
            )

    if bond_order_changes:
        print(f"\nBond Order Changes ({len(bond_order_changes)} bonds):")
        for (i, j), (order1, order2) in bond_order_changes.items():
            s1 = g1.nodes[i].get("symbol", "?")
            s2 = g1.nodes[j].get("symbol", "?")
            order1_str = interpret_bond_order(order1)
            order2_str = interpret_bond_order(order2)
            print(
                f"  Bond ({i}, {j}) [{s1}-{s2}]: "
                f"{order1_str}→{order2_str} (order {order1:.1f}→{order2:.1f})"
            )

    if charge_shift:
        print(f"\nFormal Charge Redistribution ({len(charge_shift)} atoms):")
        for i, dq in charge_shift.items():
            sym = g1.nodes[i].get("symbol", "?")
            q1 = g1.nodes[i].get("formal_charge", 0)
            q2 = g2.nodes[i].get("formal_charge", 0)
            print(f"  Atom {i} [{sym}]: charge {q1:+.0f}→{q2:+.0f} (Δq = {dq:+.2f})")

    # ASCII visualization (if generated)
    if "ascii_ts" in results:
        print("\n" + "=" * 80)
        print("ASCII REPRESENTATIONS")
        print("=" * 80)
        print("\nTransition State (TS):\n")
        print(results["ascii_ts"])
        print("\nFrame 1:\n")
        print(results.get("ascii_ref", "<no ascii_ref>"))
        print("\nFrame 2:\n")
        print(results.get("ascii_disp", "<no ascii_disp>"))
    elif debug:
        logger.debug(
            "No ASCII data available. Run with -d or ensure "
            "generate_ascii_summary() is called."
        )
        logger.debug(
            "No ASCII data available. Run with -d or ensure "
            "generate_ascii_summary() is called."
        )
