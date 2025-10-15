"""
Lightweight vibrational graph comparison using xyzgraph.

This module handles graph-based analysis of vibrational trajectories using
xyzgraph for molecular graph construction and NetworkX for graph operations.
Simplified to use distance-based detection and xyzgraph charge handling.
"""

import numpy as np
from typing import List, Tuple, Dict, Any, Optional
from ase import Atoms
from xyzgraph import build_graph, graph_to_ascii, DATA
import networkx as nx
import logging

from . import config

logger = logging.getLogger("vib_analysis")

# =====================================================================================
# === GEOMETRY HELPERS ===
# =====================================================================================

def _atoms_to_xyz_format(frame: Atoms) -> List[Tuple[str, Tuple[float, float, float]]]:
    """Convert ASE Atoms object to xyzgraph format: [(symbol, (x, y, z)), ...]"""
    return [(atom.symbol, tuple(atom.position)) for atom in frame]


def _get_distance_changes(frame_ts: Atoms, frames_displaced: List[Atoms],
                          bond: Tuple[int, int]) -> Dict[str, float]:
    """Return direct distances for TS, displaced-1, displaced-2."""
    i, j = bond
    d_ts = frame_ts.get_distance(i, j)
    d_f1 = frames_displaced[0].get_distance(i, j)
    d_f2 = frames_displaced[1].get_distance(i, j)
    return {"ts": d_ts, "f1": d_f1, "f2": d_f2}


# =====================================================================================
# === CORE GRAPH BUILDERS ===
# =====================================================================================

def build_ts_graph(
    frame_ts: Atoms,
    vib_bonds: List[Tuple[int, int]],
    vib_bond_info: Dict[Tuple[int, int], Tuple[float, float]],
    frames_displaced: List[Atoms],
    distance_tolerance: float = config.DISTANCE_TOLERANCE,
    method: str = "cheminf",
    charge: int = 0,
    multiplicity: Optional[int] = None
) -> nx.Graph:
    """
    Build transition-state reference graph using xyzgraph.
    
    Bond formation/breaking decided from direct distance changes between
    TS and displaced frames.
    
    Args:
        frame_ts: Transition state frame
        vib_bonds: List of bond tuples with significant changes
        vib_bond_info: Dict mapping bonds to (change, initial_length)
        frames_displaced: List of displaced frames (typically 2)
        distance_tolerance: Tolerance for bond formation/breaking (Å)
        method: Graph building method ('cheminf' or 'xtb')
        charge: Molecular charge
        multiplicity: Spin multiplicity (auto-detected if None)
        
    Returns:
        NetworkX graph of transition state with bond annotations
    """
    atoms_xyz = _atoms_to_xyz_format(frame_ts)
    bonds_to_add, bonds_to_remove = [], []

    for bond in vib_bonds:
        distances = _get_distance_changes(frame_ts, frames_displaced, bond)
        d_ts, d_f1, d_f2 = distances["ts"], distances["f1"], distances["f2"]

        delta, init_len = vib_bond_info.get(bond, vib_bond_info.get((bond[1], bond[0]), (0.0, 0.0)))
        avg_disp_dist = (d_f1 + d_f2) / 2

        # If both displaced distances shorter than TS by threshold → formation
        if d_ts - avg_disp_dist > distance_tolerance:
            bonds_to_add.append(bond)
        # If both displaced distances longer than TS by threshold → breaking
        elif avg_disp_dist - d_ts > distance_tolerance:
            bonds_to_remove.append(bond)

    # Build TS graph using xyzgraph; supply bond/unbond lists
    g_ts = build_graph(
        atoms_xyz,
        method=method,
        charge=charge,
        multiplicity=multiplicity,
        bond=sorted(set(bonds_to_add)),
        unbond=sorted(set(bonds_to_remove))
    )

    # Annotate vibrational bonds
    for (i, j), (delta, init_len) in vib_bond_info.items():
        if g_ts.has_edge(i, j):
            g_ts[i][j]["vib_delta"] = delta
            g_ts[i][j]["vib_initial"] = init_len
            g_ts[i][j]["vib_identified"] = True
            g_ts[i][j]["TS"] = True
        else:
            g_ts.add_edge(i, j, vib_delta=delta, vib_initial=init_len, vib_identified=True, TS=True)

    logger.debug(
        f"TS graph built: +{len(bonds_to_add)} formed, -{len(bonds_to_remove)} broken bonds."
    )
    return g_ts


def build_displaced_graphs(
    frame_ts: Atoms,
    vib_bonds: List[Tuple[int, int]],
    frames_displaced: List[Atoms],
    method: str = "cheminf",
    charge: int = 0,
    multiplicity: Optional[int] = None
) -> Tuple[nx.Graph, nx.Graph]:
    """
    Build both displaced graphs from TS geometry with guided bonding.
    
    Bonding decisions based on distance comparison between displaced frames:
    - If vib bond shorter in frame 1: add to g1's bond list, add to g2's unbond list
    - If vib bond shorter in frame 2: add to g2's bond list, add to g1's unbond list
    
    Returns:
        (g1, g2): Tuple of graphs for frames_displaced[0] and frames_displaced[1]
    """
    atoms_xyz = _atoms_to_xyz_format(frame_ts)
    
    bonds_g1 = []
    unbonds_g1 = []
    bonds_g2 = []
    unbonds_g2 = []
    
    for (i, j) in vib_bonds:
        d1 = frames_displaced[0].get_distance(i, j)
        d2 = frames_displaced[1].get_distance(i, j)
        
        if d1 < d2:  # Shorter in frame 1 → forming in g1, breaking in g2
            bonds_g1.append((i, j))
            unbonds_g2.append((i, j))
        elif d2 < d1:  # Shorter in frame 2 → breaking in g1, forming in g2
            unbonds_g1.append((i, j))
            bonds_g2.append((i, j))
    
    # Build both graphs from TS geometry with guided bonding
    g1 = build_graph(
        atoms_xyz,
        method=method,
        charge=charge,
        multiplicity=multiplicity,
        bond=sorted(set(bonds_g1)),
        unbond=sorted(set(unbonds_g1))
    )
    
    g2 = build_graph(
        atoms_xyz,
        method=method,
        charge=charge,
        multiplicity=multiplicity,
        bond=sorted(set(bonds_g2)),
        unbond=sorted(set(unbonds_g2))
    )
    
    logger.debug(f"Displaced graphs built: g1 +{len(bonds_g1)}/-{len(unbonds_g1)}, g2 +{len(bonds_g2)}/-{len(unbonds_g2)} bonds.")
    
    return g1, g2


# =====================================================================================
# === GRAPH COMPARISON AND CHARGE REDISTRIBUTION ===
# =====================================================================================

def compare_graphs(g1: nx.Graph, g2: nx.Graph) -> Dict[str, Any]:
    """
    Compare two xyzgraph graphs for bond and charge differences.
    
    Args:
        g1: First molecular graph
        g2: Second molecular graph
        
    Returns:
        Dictionary containing:
            - bonds_formed: Bonds present in g2 but not g1
            - bonds_broken: Bonds present in g1 but not g2
            - bonds_common: Bonds present in both graphs
            - bond_order_changes: Bond order changes for common bonds
            - charge_redistribution: Formal charge changes per atom
    """
    edges1, edges2 = set(g1.edges()), set(g2.edges())

    bonds_formed = list(edges2 - edges1)
    bonds_broken = list(edges1 - edges2)
    bonds_common = edges1 & edges2

    # Bond order changes for common bonds
    bond_order_changes = {}
    for (i, j) in bonds_common:
        order1 = g1[i][j].get('bond_order', 1)
        order2 = g2[i][j].get('bond_order', 1)
        if abs(order1 - order2) > 0.1:
            bond_order_changes[(i, j)] = (order1, order2)

    # Formal charge redistribution
    charge_deltas = {}
    for n in g1.nodes:
        q1 = g1.nodes[n].get("formal_charge", 0)
        q2 = g2.nodes[n].get("formal_charge", 0)
        dq = q2 - q1
        if abs(dq) > config.BOND_ORDER_EPSILON:
            charge_deltas[n] = dq

    return {
        "bonds_formed": bonds_formed,
        "bonds_broken": bonds_broken,
        "bonds_common": list(bonds_common),
        "bond_order_changes": bond_order_changes,
        "charge_redistribution": charge_deltas,
    }


# =====================================================================================
# === HIGH-LEVEL ANALYSIS PIPELINE ===
# =====================================================================================

def analyze_displacement_graphs(
        frames: List[Atoms],
        internal_changes: Dict[str, Any],
        atoms_of_interest: Optional[List[int]] = None,
        method: str = "cheminf",
        charge: int = 0,
        multiplicity: Optional[int] = None,
        distance_tolerance: float = 0.2,
        ascii_neighbor_shells: int = 1,
        ascii_scale: float = 3.0,
        ascii_include_h: bool = True,
        debug: bool = False
    ) -> Dict[str, Any]:
    """
    Full lightweight analysis: build TS + displaced graphs, classify bond and charge changes.
    
    Args:
        atoms_of_interest: Optional list of atom indices to highlight (for rotations/inversions).
                          If provided, these atoms will be used as the core for ASCII visualization
                          instead of reactive_atoms from bond changes.
    """
    ts_idx = internal_changes.get("ts_frame", 0)
    f1_idx, f2_idx = internal_changes["frame_indices"]
    vib_bond_info = internal_changes["bond_changes"]
    vib_bonds = list(vib_bond_info.keys())

    frame_ts = frames[ts_idx]
    frames_disp = [frames[f1_idx], frames[f2_idx]]

    g_ts = build_ts_graph(frame_ts, vib_bonds, vib_bond_info,
                          frames_disp, distance_tolerance,
                          method=method, charge=charge, multiplicity=multiplicity)

    # Build both displaced graphs from TS geometry with guided bonding
    g1, g2 = build_displaced_graphs(frame_ts, vib_bonds, frames_disp,
                                     method=method, charge=charge, multiplicity=multiplicity)

    comparison = compare_graphs(g1, g2)

    # Use atoms_of_interest if provided (for rotations/inversions),
    # otherwise use reactive_atoms from bond changes
    is_rotation_or_inversion = False
    if atoms_of_interest:
        reactive_atoms = set(atoms_of_interest)
        is_rotation_or_inversion = True
        
        # Mark edges between atoms_of_interest with "TS" property for highlighting
        for i in atoms_of_interest:
            for j in atoms_of_interest:
                if i < j and g_ts.has_edge(i, j):
                    g_ts[i][j]["TS"] = True
    else:
        reactive_atoms = set()
        for (i, j) in vib_bonds:
            reactive_atoms.update([i, j])

    # Expand by N shells (use union of neighbors from all graphs to stay consistent)
    for _ in range(ascii_neighbor_shells):
        all_neighbors = set()
        for n in reactive_atoms:
            if n in g_ts:
                all_neighbors.update(g_ts.neighbors(n))
            if n in g1:
                all_neighbors.update(g1.neighbors(n))
            if n in g2:
                all_neighbors.update(g2.neighbors(n))
        reactive_atoms.update(all_neighbors)
    
    # Extract consistent subgraphs
    sub_ts = g_ts.subgraph(sorted(reactive_atoms & set(g_ts.nodes))).copy()
    sub_1 = g1.subgraph(sorted(reactive_atoms & set(g1.nodes))).copy()
    sub_2 = g2.subgraph(sorted(reactive_atoms & set(g2.nodes))).copy()

    # Use subgraphs for ASCII visualization
    # For rotations/inversions, only show TS (bonding doesn't change)
    ascii_data = generate_ascii_summary(
        sub_ts, sub_1, sub_2, 
        scale=ascii_scale, 
        include_h=ascii_include_h,
        only_ts=is_rotation_or_inversion
    )
    
    results = {
        "ts_graph": g_ts,
        "frame1_graph": g1,
        "frame2_graph": g2,
        "ts_subgraph": sub_ts,
        "frame1_subgraph": sub_1,
        "frame2_subgraph": sub_2,
        "comparison": comparison,
        **ascii_data
    }


    if debug:
        logger.debug(f"Bonds formed: {comparison['bonds_formed']}")
        logger.debug(f"Bonds broken: {comparison['bonds_broken']}")
        logger.debug(f"Charge shifts: {comparison['charge_redistribution']}")

    return results


# =====================================================================================
# === OPTIONAL ASCII VISUALIZATION ===
# =====================================================================================

def generate_ascii_summary(graph_ts: nx.Graph, graph_1: nx.Graph, graph_2: nx.Graph,
                           scale: float = 3.0, include_h: bool = True, 
                           only_ts: bool = False) -> Dict[str, str]:
    """
    Generate ASCII visualization for quick debugging.
    
    Args:
        only_ts: If True, only generate TS ASCII (for rotations/inversions where bonding doesn't change)
    """
    try:
        ascii_ts = graph_to_ascii(graph_ts, scale=scale, include_h=include_h)
        
        if only_ts:
            # For rotations/inversions, skip Frame 1/2 (bonding is same)
            return {"ascii_ts": ascii_ts}
        else:
            # For bond changes, show all three
            ascii_1 = graph_to_ascii(graph_1, scale=scale, include_h=include_h, reference=graph_ts)
            ascii_2 = graph_to_ascii(graph_2, scale=scale, include_h=include_h, reference=graph_ts)
            return {"ascii_ts": ascii_ts, "ascii_ref": ascii_1, "ascii_disp": ascii_2}
    except Exception as e:
        logger.warning(f"ASCII generation failed: {e}")
        if only_ts:
            return {"ascii_ts": "<ascii_error>"}
        else:
            return {"ascii_ts": "<ascii_error>", "ascii_ref": "<ascii_error>", "ascii_disp": "<ascii_error>"}
