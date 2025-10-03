"""Graph comparison using EXACT backends.py implementation"""
import networkx as nx
from typing import Dict, List, Tuple, Any, Optional
from ase import Atoms
from ase.data import covalent_radii, atomic_numbers
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

from .data_loader import load_vdw_radii, load_expected_valences

import logging
logger = logging.getLogger('vib_analysis')


def build_graph_from_frame(frame: Atoms, use_cheminf: bool = True) -> nx.Graph:
    """
    Build molecular graph - EXACT copy from your working backends.py
    """
    if not use_cheminf:
        return _build_graph_simple(frame)
    
    # Use EXACT backends.py compute_cheminf
    vdw_data = load_vdw_radii()
    
    # Stage 1: VDW bond detection
    bonds, initial_bond_orders = _compute_cheminf_basic(frame, vdw_data)
    
    if not bonds:
        logger.warning("No bonds detected")
        return _empty_graph(frame)
    
    # Stage 2: Aromatic rings
    aromatic_bonds = _check_aromatic_rings(bonds, initial_bond_orders, frame)
    working_bond_orders = _apply_aromatic_corrections(bonds, initial_bond_orders, aromatic_bonds)
    
    # Stage 3: Valence + distance validation (EXACT from backends.py)
    final_bond_orders = _validate_bond_orders_with_valence_distance(frame, bonds, working_bond_orders, vdw_data)
    
    # Stage 4: Charges
    charges = _compute_charges_from_final_structure(frame, bonds, final_bond_orders)
    
    # Build graph
    return _build_nx_graph(frame, bonds, final_bond_orders, charges)


def _empty_graph(frame: Atoms) -> nx.Graph:
    """Empty graph with just nodes"""
    G = nx.Graph()
    for i, atom in enumerate(frame):
        G.add_node(i, symbol=atom.symbol, atomic_number=atom.number, charges={}, agg_charge=0.0)
    return G


def _build_nx_graph(frame: Atoms, bonds: List[Tuple[int, int]], bond_orders: List[float], charges: List[float]) -> nx.Graph:
    """Construct NetworkX graph from computed data"""
    G = nx.Graph()
    
    for i, atom in enumerate(frame):
        G.add_node(i, 
                  symbol=atom.symbol, 
                  atomic_number=atom.number,
                  charges={'gasteiger': charges[i]},
                  agg_charge=charges[i])
    
    for idx, (i, j) in enumerate(bonds):
        G.add_edge(i, j,
                  bond_order=bond_orders[idx],
                  distance=frame.get_distance(i, j),
                  bond_type=(frame[i].symbol, frame[j].symbol))
    
    return G


def _build_graph_simple(frame: Atoms) -> nx.Graph:
    """Fallback simple construction"""
    G = nx.Graph()
    
    for i, atom in enumerate(frame):
        G.add_node(i, symbol=atom.symbol, atomic_number=atom.number)
    
    for i in range(len(frame)):
        for j in range(i+1, len(frame)):
            distance = frame.get_distance(i, j)
            r1 = covalent_radii[frame[i].number]
            r2 = covalent_radii[frame[j].number]
            
            if distance < 1.4 * (r1 + r2):
                G.add_edge(i, j, bond_order=1.0, distance=distance)
    
    return G


def _compute_cheminf_basic(frame: Atoms, vdw_data: Dict[str, float]) -> Tuple[List[Tuple[int, int]], List[float]]:
    """Stage 1: EXACT from backends.py"""
    bonds = []
    bond_orders = []
    
    for i in range(len(frame)):
        for j in range(i+1, len(frame)):
            distance = np.linalg.norm(frame.positions[i] - frame.positions[j])
            r1 = vdw_data.get(frame[i].symbol, 2.0)
            r2 = vdw_data.get(frame[j].symbol, 2.0)
            vdw_sum = r1 + r2
            
            # Thresholds from backends.py
            if 'H' in (frame[i].symbol, frame[j].symbol):
                threshold = 0.45 * vdw_sum
            else:
                threshold = 0.55 * vdw_sum
            
            if distance < threshold:
                # Metal check
                if _should_bond_metal(i, j, distance, frame, vdw_data):
                    bonds.append((i, j))
                    bond_orders.append(1.0)  # Always start at 1.0
    
    return bonds, bond_orders


def _should_bond_metal(i: int, j: int, distance: float, frame: Atoms, vdw_data: Dict[str, float]) -> bool:
    """EXACT from backends.py"""
    metals = ('Li', 'Na', 'K', 'Mg', 'Ca', 'Zn', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
              'Ru', 'Rh', 'Pd', 'Ag', 'Ir', 'Pt', 'Au')
    
    sym_i, sym_j = frame[i].symbol, frame[j].symbol
    
    if sym_i not in metals and sym_j not in metals:
        return True
    
    metal_sym = sym_i if sym_i in metals else sym_j
    other_sym = sym_j if sym_i in metals else sym_i
    
    r_metal = vdw_data.get(metal_sym, 2.0)
    r_other = vdw_data.get(other_sym, 2.0)
    vdw_contact = r_metal + r_other
    
    # Ionic bonding
    if ((metal_sym in ('Li', 'Na', 'K') and other_sym in ('F', 'Cl', 'Br', 'I')) or
        (metal_sym in ('Mg', 'Ca', 'Zn') and other_sym in ('F', 'Cl', 'Br', 'I', 'O')) or
        (metal_sym in ('Mg', 'Ca') and other_sym == 'S')):
        return distance < 0.4 * vdw_contact
    
    # Coordination
    if other_sym in ('O', 'N', 'C', 'P', 'S') and distance < 0.4 * vdw_contact:
        return True
    
    return False


def _check_aromatic_rings(bonds: List[Tuple[int, int]], bond_orders: List[float], frame: Atoms) -> List[Tuple[int, int]]:
    """EXACT from backends.py"""
    graph = nx.Graph()
    graph.add_nodes_from(range(len(frame)))
    for i, j in bonds:
        graph.add_edge(i, j)
    
    cycles = nx.cycle_basis(graph)
    aromatic_bonds = []
    
    for cycle in cycles:
        if len(cycle) in [5, 6]:
            if all(frame[atom].symbol in ['C', 'N'] for atom in cycle):
                for k in range(len(cycle)):
                    atom1 = cycle[k]
                    atom2 = cycle[(k + 1) % len(cycle)]
                    bond_pair = (min(atom1, atom2), max(atom1, atom2))
                    if bond_pair in [(i, j) for i, j in bonds]:
                        aromatic_bonds.append(bond_pair)
    
    return aromatic_bonds


def _apply_aromatic_corrections(bonds: List[Tuple[int, int]], bond_orders: List[float], aromatic_bonds: List[Tuple[int, int]]) -> List[float]:
    """EXACT from backends.py"""
    final_orders = bond_orders.copy()
    
    for idx, (i, j) in enumerate(bonds):
        bond_pair = (min(i, j), max(i, j))
        if bond_pair in aromatic_bonds:
            final_orders[idx] = 1.5
    
    return final_orders


def _validate_bond_orders_with_valence_distance(frame: Atoms, bonds: List[Tuple[int, int]], 
                                                initial_orders: List[float], vdw_data: Dict[str, float]) -> List[float]:
    """EXACT from backends.py validate_bond_orders_with_valence_distance"""
    expected_valences = load_expected_valences()
    final_orders = initial_orders.copy()
    
    max_iterations = 5
    for iteration in range(max_iterations):
        changes = 0
        
        current_graph = _build_temp_graph(frame, bonds, final_orders)
        problematic_atoms = _identify_valence_problems(current_graph, expected_valences)
        
        if not problematic_atoms:
            break
        
        for atom_idx, valence_info in problematic_atoms.items():
            atom_changes = _fix_atom_valence(frame, bonds, final_orders, atom_idx, 
                                            valence_info, problematic_atoms, vdw_data, expected_valences)
            changes += atom_changes
        
        if changes == 0:
            break
    
    return final_orders


def _build_temp_graph(frame: Atoms, bonds: List[Tuple[int, int]], bond_orders: List[float]) -> nx.Graph:
    """Build temporary graph"""
    G = nx.Graph()
    for i, atom in enumerate(frame):
        G.add_node(i, symbol=atom.symbol, atomic_number=atom.number)
    
    for idx, (i, j) in enumerate(bonds):
        G.add_edge(i, j, bond_order=bond_orders[idx])
    
    return G


def _identify_valence_problems(graph: nx.Graph, expected_valences: Dict[str, List[int]]) -> Dict[int, Dict[str, Any]]:
    """EXACT from backends.py"""
    problematic = {}
    
    for i, node_data in graph.nodes(data=True):
        sym = node_data.get('symbol', '')
        
        if sym not in expected_valences:
            continue
        
        total_bond_order = sum(
            graph.edges[i, neighbor].get('bond_order', 1.0) 
            for neighbor in graph.neighbors(i)
        )
        
        expected = expected_valences[sym]
        valence_error = min(abs(total_bond_order - exp) for exp in expected)
        
        if valence_error > 0.5:
            problematic[i] = {
                'current_valence': total_bond_order,
                'expected_valences': expected,
                'error': valence_error,
                'symbol': sym
            }
    
    return problematic


def _fix_atom_valence(frame: Atoms, bonds: List[Tuple[int, int]], bond_orders: List[float],
                      atom_idx: int, valence_info: Dict[str, Any], 
                      all_problematic: Dict[int, Dict[str, Any]], vdw_data: Dict[str, float],
                      expected_valences: Dict[str, List[int]]) -> int:
    """EXACT from backends.py fix_atom_valence"""
    changes = 0
    current_valence = valence_info['current_valence']
    expected = valence_info['expected_valences']
    
    target_valence = min(expected, key=lambda x: abs(x - current_valence))
    valence_deficit = target_valence - current_valence
    
    if abs(valence_deficit) < 0.3:
        return 0
    
    candidate_bonds = []
    for bond_idx, (i, j) in enumerate(bonds):
        if i == atom_idx or j == atom_idx:
            other_atom = j if i == atom_idx else i
            
            if 'H' in (frame[atom_idx].symbol, frame[other_atom].symbol):
                continue
            
            if other_atom not in all_problematic:
                continue
            
            distance = np.linalg.norm(frame.positions[atom_idx] - frame.positions[other_atom])
            r1 = vdw_data.get(frame[atom_idx].symbol, 2.0)
            r2 = vdw_data.get(frame[other_atom].symbol, 2.0)
            ratio = distance / (r1 + r2)
            current_order = bond_orders[bond_idx]
            
            other_valence_info = all_problematic[other_atom]
            other_target = min(other_valence_info['expected_valences'], 
                             key=lambda x: abs(x - other_valence_info['current_valence']))
            other_deficit = other_target - other_valence_info['current_valence']
            
            # THREE CASES from backends.py
            is_complementary = False
            potential_increase = 0
            
            if valence_deficit > 0.3 and other_deficit > 0.3:
                potential_increase = min(valence_deficit, other_deficit)
                is_complementary = True
            elif valence_deficit > 0.3 and other_deficit < -0.3:
                potential_increase = min(valence_deficit, abs(other_deficit))
                is_complementary = True
            elif valence_deficit < -0.3 and other_deficit > 0.3:
                potential_increase = min(abs(valence_deficit), other_deficit)
                is_complementary = True
            
            if is_complementary and ratio < 0.55:
                max_order = 3.0 if ratio < 0.35 else 2.0
                final_increase = min(potential_increase, max_order - current_order)
                
                if final_increase >= 0.5:
                    candidate_bonds.append((bond_idx, other_atom, final_increase, ratio))
    
    if not candidate_bonds:
        return 0
    
    candidate_bonds.sort(key=lambda x: x[2], reverse=True)
    bond_idx, other_atom, increase, ratio = candidate_bonds[0]
    
    bond_orders[bond_idx] += increase
    changes += 1
    
    return changes


def _compute_charges_from_final_structure(frame: Atoms, bonds: List[Tuple[int, int]], bond_orders: List[float]) -> List[float]:
    """EXACT from backends.py"""
    try:
        mol = Chem.RWMol()
        
        for i, atom in enumerate(frame):
            rd_atom = Chem.Atom(atom.symbol)
            mol.AddAtom(rd_atom)
        
        for idx, (i, j) in enumerate(bonds):
            bo = bond_orders[idx]
            if bo >= 2.5:
                bt = Chem.BondType.TRIPLE
            elif bo >= 1.75:
                bt = Chem.BondType.DOUBLE
            elif bo >= 1.25:
                bt = Chem.BondType.AROMATIC
            else:
                bt = Chem.BondType.SINGLE
            
            mol.AddBond(int(i), int(j), bt)
        
        Chem.SanitizeMol(mol)
        AllChem.ComputeGasteigerCharges(mol)
        
        charges = []
        for atom in mol.GetAtoms():
            charge = atom.GetProp('_GasteigerCharge')
            if charge == 'nan' or (isinstance(charge, float) and np.isnan(charge)):
                charges.append(0.0)
            else:
                charges.append(float(charge))
        
        return charges
        
    except Exception as e:
        logger.warning(f"Charge computation failed: {e}")
        return [0.0] * len(frame)


def validate_graph_valences(graph: nx.Graph, frame: Atoms, expected_valences: Dict[str, List[int]]) -> Dict[str, Any]:
    """Validate graph valences and report issues"""
    report = {
        'total_atoms': graph.number_of_nodes(),
        'total_bonds': graph.number_of_edges(),
        'valence_errors': {},
        'hypervalent_atoms': [],
        'formal_charges': {},
        'bond_order_summary': {},
        'charge_distribution': {}
    }
    
    # Collect bond order statistics
    for i, j, data in graph.edges(data=True):
        bo = data.get('bond_order', 1.0)
        report['bond_order_summary'][bo] = report['bond_order_summary'].get(bo, 0) + 1
    
    # Check each atom's valence and charges
    for i, node_data in graph.nodes(data=True):
        sym = node_data.get('symbol', '')
        formal_charge = node_data.get('formal_charge', 0)
        gasteiger = node_data.get('charges', {}).get('gasteiger', 0.0)
        
        report['charge_distribution'][i] = {
            'symbol': sym,
            'formal_charge': formal_charge,
            'gasteiger': gasteiger
        }
        
        if formal_charge != 0:
            report['formal_charges'][i] = formal_charge
        
        if sym not in expected_valences:
            continue
        
        total_bond_order = sum(
            graph.edges[i, neighbor].get('bond_order', 1.0) 
            for neighbor in graph.neighbors(i)
        )
        
        expected = expected_valences[sym]
        valence_error = min(abs(total_bond_order - exp) for exp in expected)
        
        if valence_error > 0.5:
            report['valence_errors'][i] = {
                'symbol': sym,
                'current_valence': total_bond_order,
                'expected': expected,
                'error': valence_error,
                'neighbors': list(graph.neighbors(i)),
                'bonds': [(i, n, graph.edges[i, n].get('bond_order', 1.0)) for n in graph.neighbors(i)]
            }
        
        max_expected = max(expected)
        if total_bond_order > max_expected + 0.5:
            report['hypervalent_atoms'].append({
                'atom_idx': i,
                'symbol': sym,
                'valence': total_bond_order,
                'max_expected': max_expected,
                'excess': total_bond_order - max_expected
            })
    
    return report


def print_graph_validation(report: Dict[str, Any], frame: Atoms):
    """Print detailed validation report"""
    print("\n" + "="*80)
    print(" "*25 + "GRAPH VALIDATION REPORT")
    print("="*80)
    
    print(f"\nTotal atoms: {report['total_atoms']}")
    print(f"Total bonds: {report['total_bonds']}")
    
    if report['bond_order_summary']:
        print("\nBond order distribution:")
        for bo, count in sorted(report['bond_order_summary'].items()):
            print(f"  {bo:.1f}: {count} bonds")
    
    if report['charge_distribution']:
        charges = [info['gasteiger'] for info in report['charge_distribution'].values()]
        print(f"\nCharge statistics:")
        print(f"  Mean: {np.mean(charges):.4f}, Std: {np.std(charges):.4f}")
        print(f"  Min: {np.min(charges):.4f}, Max: {np.max(charges):.4f}")
    
    if report['valence_errors']:
        print(f"\n⚠ Valence errors: {len(report['valence_errors'])} atoms")
        for idx, info in report['valence_errors'].items():
            print(f"  Atom {idx} ({info['symbol']}): valence={info['current_valence']:.2f}, expected={info['expected']}")
    
    if report['hypervalent_atoms']:
        print(f"\n⚠ Hypervalent atoms: {len(report['hypervalent_atoms'])}")
        for hv in report['hypervalent_atoms']:
            print(f"  Atom {hv['atom_idx']} ({hv['symbol']}): {hv['valence']:.2f} (max={hv['max_expected']})")
    
    if not report['valence_errors'] and not report['hypervalent_atoms']:
        print("\n✓ All valences within expected ranges")
    
    print("="*80)


def cross_validate_with_internal_coords(graph: nx.Graph, internal_bond_changes: Dict[Tuple[int, int], Tuple[float, float]], 
                                       frame: Atoms) -> Dict[str, Any]:
    """Cross-validate graph with vibrational bond changes"""
    graph_bonds = set(graph.edges())
    changing_bonds = set(internal_bond_changes.keys())
    
    validation = {
        'graph_detects_changes': [],
        'graph_misses_changes': []
    }
    
    for bond in changing_bonds:
        norm_bond = tuple(sorted(bond))
        if norm_bond in graph_bonds or tuple(reversed(norm_bond)) in graph_bonds:
            change_mag, initial = internal_bond_changes[bond]
            bo = graph.edges.get(norm_bond, {}).get('bond_order', 
                                 graph.edges.get(tuple(reversed(norm_bond)), {}).get('bond_order', 1.0))
            validation['graph_detects_changes'].append({
                'bond': bond,
                'change': change_mag,
                'initial': initial,
                'bond_order': bo
            })
        else:
            validation['graph_misses_changes'].append({
                'bond': bond,
                'change': internal_bond_changes[bond][0],
                'distance': frame.get_distance(*bond)
            })
    
    return validation


def compare_graphs(ref_graph: nx.Graph, displaced_graph: nx.Graph, 
                   frame_ref: Atoms, frame_disp: Atoms) -> Dict[str, Any]:
    """Compare two molecular graphs"""
    ref_edges = set(ref_graph.edges())
    disp_edges = set(displaced_graph.edges())
    
    bonds_formed = list(disp_edges - ref_edges)
    bonds_broken = list(ref_edges - disp_edges)
    
    # Bond order changes
    common_bonds = ref_edges & disp_edges
    bond_order_changes = {}
    
    for i, j in common_bonds:
        ref_order = ref_graph[i][j].get('bond_order', 1.0)
        disp_order = displaced_graph[i][j].get('bond_order', 1.0)
        
        if abs(ref_order - disp_order) >= 0.3:
            bond_order_changes[(i, j)] = (ref_order, disp_order)
    
    # Ring changes
    ref_rings = find_rings(ref_graph)
    disp_rings = find_rings(displaced_graph)
    
    rings_formed = [r for r in disp_rings if r not in ref_rings]
    rings_broken = [r for r in ref_rings if r not in disp_rings]
    
    # Fragments
    ref_fragments = nx.number_connected_components(ref_graph)
    disp_fragments = nx.number_connected_components(displaced_graph)
    
    # Classify
    transformation_type = classify_transformation(
        bonds_formed, bonds_broken, bond_order_changes,
        rings_formed, rings_broken, ref_fragments, disp_fragments
    )
    
    return {
        'bonds_formed': bonds_formed,
        'bonds_broken': bonds_broken,
        'bond_order_changes': bond_order_changes,
        'rings_formed': rings_formed,
        'rings_broken': rings_broken,
        'ref_fragments': ref_fragments,
        'disp_fragments': disp_fragments,
        'transformation_type': transformation_type,
    }


def find_rings(graph: nx.Graph, max_size: int = 12) -> List[Tuple[int, ...]]:
    """Find all rings up to max_size"""
    try:
        cycles = nx.cycle_basis(graph)
        rings = [tuple(sorted(cycle)) for cycle in cycles if len(cycle) <= max_size]
        return rings
    except:
        return []


def classify_transformation(bonds_formed, bonds_broken, bond_order_changes,
                           rings_formed, rings_broken, ref_frags, disp_frags) -> str:
    """Classify transformation type"""
    classifications = []
    
    if disp_frags > ref_frags:
        classifications.append("fragmentation")
    elif disp_frags < ref_frags:
        classifications.append("bond_formation")
    
    if rings_formed and not rings_broken:
        classifications.append("cyclization")
    elif rings_broken and not rings_formed:
        classifications.append("ring_opening")
    elif rings_formed and rings_broken:
        classifications.append("ring_rearrangement")
    
    if bond_order_changes:
        increasing = sum(1 for old, new in bond_order_changes.values() if new > old)
        decreasing = sum(1 for old, new in bond_order_changes.values() if new < old)
        
        if increasing and decreasing:
            classifications.append("bond_order_redistribution")
        elif increasing:
            classifications.append("bond_strengthening")
        elif decreasing:
            classifications.append("bond_weakening")
    
    if len(bonds_formed) >= 2 and len(bonds_broken) >= 2:
        classifications.append("concerted_rearrangement")
    
    if (bonds_formed or bonds_broken) and ref_frags == disp_frags and not classifications:
        classifications.append("rearrangement")
    
    return " + ".join(classifications) if classifications else "no_significant_change"


def print_graph_analysis(graph_results: Dict[str, Any], atom_index_map: Dict[int, str] = None):
    """Print formatted graph comparison results"""
    print("\n" + "="*80)
    print(" "*25 + "GRAPH ANALYSIS SUMMARY")
    print("="*80)
    
    f1_stats = graph_results['frame1_stats']
    f2_stats = graph_results['frame2_stats']
    
    print(f"\nFrame {graph_results['frame1_idx']}:")
    print(f"  Atoms: {f1_stats['n_atoms']}, Bonds: {f1_stats['n_bonds']}")
    print(f"  Rings: {f1_stats['n_rings']}, Fragments: {f1_stats['n_fragments']}")
    if f1_stats.get('vib_augmented_bonds', 0) > 0:
        print(f"  Vib-augmented bonds: {f1_stats['vib_augmented_bonds']}")
    
    print(f"\nFrame {graph_results['frame2_idx']}:")
    print(f"  Atoms: {f2_stats['n_atoms']}, Bonds: {f2_stats['n_bonds']}")
    print(f"  Rings: {f2_stats['n_rings']}, Fragments: {f2_stats['n_fragments']}")
    if f2_stats.get('vib_augmented_bonds', 0) > 0:
        print(f"  Vib-augmented bonds: {f2_stats['vib_augmented_bonds']}")
    
    comp = graph_results['comparison']
    
    print(f"\n{'─'*80}")
    print(f"Transformation: {comp['transformation_type']}")
    
    if comp['bonds_formed']:
        print(f"\nBonds Formed ({len(comp['bonds_formed'])}):")
        for (i, j), (sym1, sym2) in zip(comp['bonds_formed'], comp['bonds_formed_symbols']):
            print(f"  ({i}, {j}) [{sym1}-{sym2}]")
    
    if comp['bonds_broken']:
        print(f"\nBonds Broken ({len(comp['bonds_broken'])}):")
        for (i, j), (sym1, sym2) in zip(comp['bonds_broken'], comp['bonds_broken_symbols']):
            print(f"  ({i}, {j}) [{sym1}-{sym2}]")
    
    if comp['bond_order_changes']:
        print(f"\nBond Order Changes ({len(comp['bond_order_changes'])}):")
        for (i, j), (old, new) in comp['bond_order_changes'].items():
            if atom_index_map:
                sym1, sym2 = atom_index_map[i], atom_index_map[j]
                print(f"  ({i}, {j}) [{sym1}-{sym2}]: {old:.1f} → {new:.1f}")
            else:
                print(f"  ({i}, {j}): {old:.1f} → {new:.1f}")
    
    if comp['rings_formed']:
        print(f"\nRings Formed: {len(comp['rings_formed'])}")
    
    if comp['rings_broken']:
        print(f"\nRings Broken: {len(comp['rings_broken'])}")
    
    print("="*80)


def augment_graph_with_missed_bonds(graph: nx.Graph, frame: Atoms, missed_bonds: List[Dict[str, Any]], 
                                     vdw_data: Dict[str, float]) -> nx.Graph:
    """
    Add bonds that were detected by vibrational analysis but missed by graph construction.
    Uses distance and VDW radii to estimate appropriate bond order.
    """
    for miss_info in missed_bonds:
        i, j = miss_info['bond']
        current_dist = miss_info['distance']
        
        sym_i, sym_j = frame[i].symbol, frame[j].symbol
        r_i = vdw_data.get(sym_i, 2.0)
        r_j = vdw_data.get(sym_j, 2.0)
        vdw_sum = r_i + r_j
        
        ratio = current_dist / vdw_sum
        
        # Conservative bond order estimation
        if ratio < 0.40:
            bond_order = 2.5
        elif ratio < 0.50:
            bond_order = 1.5
        elif ratio < 0.65:
            bond_order = 1.0
        else:
            bond_order = 0.5
        
        graph.add_edge(i, j, 
                      bond_order=bond_order, 
                      distance=current_dist,
                      bond_type=(sym_i, sym_j),
                      vib_augmented=True)
        
        logger.debug(f"Augmented: {i}-{j} ({sym_i}-{sym_j}), order={bond_order:.2f}, dist={current_dist:.2f}Å")
    
    return graph


def analyze_displacement_graphs(frames: List[Atoms], frame_indices: List[int], 
                                bond_tolerance: float = 1.4, debug: bool = False,
                                internal_changes: Dict[str, Any] = None) -> Dict[str, Any]:
    """
    Compare graphs with intelligent bond augmentation from vibrational analysis.
    """
    expected_valences = load_expected_valences()
    vdw_data = load_vdw_radii()
    
    if len(frame_indices) < 2:
        logger.warning("Need at least 2 frames")
        return {}
    
    frame1_idx, frame2_idx = frame_indices[0], frame_indices[1]
    frame1 = frames[frame1_idx]
    frame2 = frames[frame2_idx]
    
    # Build initial graphs
    graph1 = build_graph_from_frame(frame1, use_cheminf=True)
    graph2 = build_graph_from_frame(frame2, use_cheminf=True)
    
    # Intelligent augmentation if vib data provided
    if internal_changes and 'bond_changes' in internal_changes:
        vib_bond_changes = internal_changes['bond_changes']
        
        if debug:
            print(f"\n{'='*80}")
            print("Vibrational Analysis Integration")
            print(f"{'='*80}")
            print(f"Detected {len(vib_bond_changes)} changing bonds from vibration")
        
        # Check which bonds are missed
        cv1 = cross_validate_with_internal_coords(graph1, vib_bond_changes, frame1)
        cv2 = cross_validate_with_internal_coords(graph2, vib_bond_changes, frame2)
        
        if debug and (cv1['graph_misses_changes'] or cv2['graph_misses_changes']):
            print(f"\nFrame {frame1_idx} missing {len(cv1['graph_misses_changes'])} bonds")
            print(f"Frame {frame2_idx} missing {len(cv2['graph_misses_changes'])} bonds")
        
        # Augment based on distance comparison
        for (i, j), (change_mag, initial_dist) in vib_bond_changes.items():
            dist1 = frame1.get_distance(i, j)
            dist2 = frame2.get_distance(i, j)
            
            norm_bond = tuple(sorted((i, j)))
            in_graph1 = graph1.has_edge(*norm_bond)
            in_graph2 = graph2.has_edge(*norm_bond)
            
            if not in_graph1 and not in_graph2:
                # Add to frame with shorter distance (more bonded state)
                if dist1 < dist2:
                    if debug:
                        print(f"Adding {i}-{j} to frame {frame1_idx} (shorter: {dist1:.2f}Å)")
                    graph1 = augment_graph_with_missed_bonds(
                        graph1, frame1, [{'bond': (i, j), 'distance': dist1}], vdw_data
                    )
                else:
                    if debug:
                        print(f"Adding {i}-{j} to frame {frame2_idx} (shorter: {dist2:.2f}Å)")
                    graph2 = augment_graph_with_missed_bonds(
                        graph2, frame2, [{'bond': (i, j), 'distance': dist2}], vdw_data
                    )
            elif not in_graph1 and dist1 < dist2 * 1.1:
                if debug:
                    print(f"Adding {i}-{j} to frame {frame1_idx} (comparable)")
                graph1 = augment_graph_with_missed_bonds(
                    graph1, frame1, [{'bond': (i, j), 'distance': dist1}], vdw_data
                )
            elif not in_graph2 and dist2 < dist1 * 1.1:
                if debug:
                    print(f"Adding {i}-{j} to frame {frame2_idx} (comparable)")
                graph2 = augment_graph_with_missed_bonds(
                    graph2, frame2, [{'bond': (i, j), 'distance': dist2}], vdw_data
                )
    
    # Debug validation
    if debug:
        print(f"\n{'='*80}")
        print(f"Post-Augmentation: Frame {frame1_idx}")
        print(f"{'='*80}")
        report1 = validate_graph_valences(graph1, frame1, expected_valences)
        print_graph_validation(report1, frame1)
        
        print(f"\n{'='*80}")
        print(f"Post-Augmentation: Frame {frame2_idx}")
        print(f"{'='*80}")
        report2 = validate_graph_valences(graph2, frame2, expected_valences)
        print_graph_validation(report2, frame2)
    
    # Compare
    comparison = compare_graphs(graph1, graph2, frame1, frame2)
    
    comparison['bonds_formed_symbols'] = [
        (frame1[i].symbol, frame1[j].symbol) for i, j in comparison['bonds_formed']
    ]
    comparison['bonds_broken_symbols'] = [
        (frame1[i].symbol, frame1[j].symbol) for i, j in comparison['bonds_broken']
    ]
    
    results = {
        'frame1_idx': frame1_idx,
        'frame2_idx': frame2_idx,
        'frame1_stats': {
            'n_atoms': len(frame1),
            'n_bonds': graph1.number_of_edges(),
            'n_rings': len(find_rings(graph1)),
            'n_fragments': nx.number_connected_components(graph1),
            'vib_augmented_bonds': len([(i,j) for i,j,d in graph1.edges(data=True) if d.get('vib_augmented')])
        },
        'frame2_stats': {
            'n_atoms': len(frame2),
            'n_bonds': graph2.number_of_edges(),
            'n_rings': len(find_rings(graph2)),
            'n_fragments': nx.number_connected_components(graph2),
            'vib_augmented_bonds': len([(i,j) for i,j,d in graph2.edges(data=True) if d.get('vib_augmented')])
        },
        'comparison': comparison
    }
    
    # Final cross-validation
    if internal_changes and 'bond_changes' in internal_changes:
        results['cross_validation'] = cross_validate_with_internal_coords(
            graph1, internal_changes['bond_changes'], frame1
        )
        
        if debug:
            print(f"\n{'='*80}")
            print("Final Cross-Validation")
            print(f"{'='*80}")
            cv = results['cross_validation']
            print(f"✓ Detects: {len(cv['graph_detects_changes'])} / {len(internal_changes['bond_changes'])} bonds")
            
            if cv['graph_misses_changes']:
                print(f"\n⚠ Still missing {len(cv['graph_misses_changes'])} bonds:")
                for miss in cv['graph_misses_changes']:
                    bond = miss['bond']
                    print(f"  {bond} ({frame1[bond[0]].symbol}-{frame1[bond[1]].symbol}): "
                          f"Δ={miss['change']:.3f}Å, d={miss['distance']:.3f}Å")
            else:
                print("\n✓ All vib-detected bonds now in graph!")
    
    return results


def print_graph_analysis(graph_results: Dict[str, Any], atom_index_map: Dict[int, str] = None):
    """Print formatted graph comparison results"""
    print("\n" + "="*80)
    print(" "*25 + "GRAPH ANALYSIS SUMMARY")
    print("="*80)
    
    f1_stats = graph_results['frame1_stats']
    f2_stats = graph_results['frame2_stats']
    
    print(f"\nFrame {graph_results['frame1_idx']}:")
    print(f"  Atoms: {f1_stats['n_atoms']}, Bonds: {f1_stats['n_bonds']}")
    print(f"  Rings: {f1_stats['n_rings']}, Fragments: {f1_stats['n_fragments']}")
    if f1_stats.get('vib_augmented_bonds', 0) > 0:
        print(f"  Vib-augmented: {f1_stats['vib_augmented_bonds']} bonds")
    
    print(f"\nFrame {graph_results['frame2_idx']}:")
    print(f"  Atoms: {f2_stats['n_atoms']}, Bonds: {f2_stats['n_bonds']}")
    print(f"  Rings: {f2_stats['n_rings']}, Fragments: {f2_stats['n_fragments']}")
    if f2_stats.get('vib_augmented_bonds', 0) > 0:
        print(f"  Vib-augmented: {f2_stats['vib_augmented_bonds']} bonds")
    
    comp = graph_results['comparison']
    
    print(f"\n{'─'*80}")
    print(f"Transformation: {comp['transformation_type']}")
    
    if comp['bonds_formed']:
        print(f"\nBonds Formed ({len(comp['bonds_formed'])}):")
        for (i, j), (sym1, sym2) in zip(comp['bonds_formed'], comp['bonds_formed_symbols']):
            print(f"  ({i}, {j}) [{sym1}-{sym2}]")
    
    if comp['bonds_broken']:
        print(f"\nBonds Broken ({len(comp['bonds_broken'])}):")
        for (i, j), (sym1, sym2) in zip(comp['bonds_broken'], comp['bonds_broken_symbols']):
            print(f"  ({i}, {j}) [{sym1}-{sym2}]")
    
    if comp['bond_order_changes']:
        print(f"\nBond Order Changes ({len(comp['bond_order_changes'])}):")
        for (i, j), (old, new) in comp['bond_order_changes'].items():
            if atom_index_map:
                sym1, sym2 = atom_index_map[i], atom_index_map[j]
                print(f"  ({i}, {j}) [{sym1}-{sym2}]: {old:.1f} → {new:.1f}")
            else:
                print(f"  ({i}, {j}): {old:.1f} → {new:.1f}")
    
    if comp.get('rings_formed'):
        print(f"\nRings Formed: {len(comp['rings_formed'])}")
    
    if comp.get('rings_broken'):
        print(f"\nRings Broken: {len(comp['rings_broken'])}")
    
    print("="*80)