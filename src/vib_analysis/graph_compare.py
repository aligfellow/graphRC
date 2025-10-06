"""Graph comparison using xyzgraph library"""
import networkx as nx
from typing import Dict, List, Tuple, Any, Optional
from ase import Atoms
import numpy as np
import itertools
import re

try:
    from xyzgraph import build_graph, graph_to_ascii, VDW as XYZ_VDW
except ImportError:
    raise ImportError("xyzgraph library required. Install with: pip install xyzgraph")

import logging
logger = logging.getLogger('vib_analysis')


# ===================== SIMPLE GEOMETRY UTILITIES =====================

def _vdw_radius(sym: str) -> float:
    return XYZ_VDW.get(sym, 2.0)

def _vdw_ratio(frame: Atoms, i: int, j: int, vdw: Dict[str,float]) -> Tuple[float,float,float]:
    dist = frame.get_distance(i, j)
    r_sum = vdw.get(frame[i].symbol, 2.0) + vdw.get(frame[j].symbol, 2.0)
    ratio = dist / r_sum if r_sum > 1e-6 else 2.0
    return ratio, dist, r_sum

def _bond_threshold(sym1: str, sym2: str) -> float:
    if 'H' in (sym1, sym2):
        return 0.45
    if sym1 in ('Ru','Mn','Fe','Co','Ni','Cu','Zn') or sym2 in ('Ru','Mn','Fe','Co','Ni','Cu','Zn'):
        return 0.70
    if sym1 in ('C','N','O') and sym2 in ('C','N','O'):
        return 0.55
    return 0.60

# ===================== TS VIB BOND SELECTION (LIGHT SCORING) =====================

def _prepare_vib_candidates(base_graph: nx.Graph,
                            frame_ts: Atoms,
                            vib_bonds: List[Tuple[int,int]],
                            vdw: Dict[str,float],
                            vib_bond_info: Dict[Tuple[int,int], Tuple[float,float]]) -> List[Dict[str,Any]]:
    # NOTE: now includes vibrational delta for scoring
    candidates = []
    for i, j in vib_bonds:
        ratio, dist, _ = _vdw_ratio(frame_ts, i, j, vdw)
        delta, _init = vib_bond_info.get((i,j), vib_bond_info.get((j,i),(0.0,0.0)))
        if base_graph.has_edge(i, j):
            candidates.append({
                "bond": (i, j),
                "already": True,
                "ratio": ratio,
                "dist": dist,
                "delta": delta
            })
        else:
            if ratio < 0.45: bo = 1.5
            elif ratio < 0.60: bo = 1.0
            elif ratio < 0.75: bo = 0.5
            else: bo = 0.3
            candidates.append({
                "bond": (i, j),
                "already": False,
                "ratio": ratio,
                "dist": dist,
                "bo": bo,
                "delta": delta
            })
    return candidates

def _score_subset(base_graph: nx.Graph,
                  frame_ts: Atoms,
                  subset: List[Dict[str,Any]],
                  vdw: Dict[str,float]) -> float:
    """
    Simple scoring:
      +10 per new 5/6 ring formed when adding subset bonds
      +3 per bond whose ratio < (threshold - 0.05) (strong geometry)
      -5 per bond with ratio > (threshold + 0.10) (weak candidate)
    """
    test_g = base_graph.copy()
    for c in subset:
        if not c["already"]:
            i, j = c["bond"]
            test_g.add_edge(i, j)
    base_rings = set(tuple(sorted(r)) for r in find_rings(base_graph))
    test_rings = set(tuple(sorted(r)) for r in find_rings(test_g))
    new_rings = test_rings - base_rings
    score = 0.0
    # ring reward (unchanged)
    if new_rings:
        for ring in new_rings:
            ln = len(ring)
            if ln in (5, 6):
                score += 10.0
    # geometric quality + delta weighting + redundancy penalty
    for c in subset:
        i, j = c["bond"]
        thr = _bond_threshold(frame_ts[i].symbol, frame_ts[j].symbol)
        ratio = c["ratio"]
        delta = abs(c.get("delta", 0.0))
        # delta normalization (simple)
        score += 5.0 * delta  # prioritize largest vibrational displacement candidates
        if ratio < thr - 0.05:
            score += 3.0
        elif ratio > thr + 0.10:
            score -= 5.0
        # redundancy penalty: very short existing path in base graph
        if i in base_graph and j in base_graph and nx.has_path(base_graph, i, j):
            try:
                sp_len = nx.shortest_path_length(base_graph, i, j)
                if sp_len <= 2:
                    score -= 8.0  # penalize “shortcut” like 13–16–15
            except Exception:
                pass
    return score

def _select_vib_subset(base_graph: nx.Graph,
                       frame_ts: Atoms,
                       candidates: List[Dict[str,Any]],
                       vdw: Dict[str,float],
                       enable_scoring: bool) -> List[Dict[str,Any]]:
    if not enable_scoring or len(candidates) == 0:
        return candidates  # take all
    # build subsets: empty, singles, pairs (if <=8), full
    subsets = [[]]
    for c in candidates:
        subsets.append([c])
    if len(candidates) <= 8 and len(candidates) >= 2:
        import itertools
        for a, b in itertools.combinations(candidates, 2):
            subsets.append([a, b])
    if len(candidates) <= 10:
        subsets.append(candidates)
    best = (None, -1e9)
    for subset in subsets:
        if not subset:
            continue
        score = _score_subset(base_graph, frame_ts, subset, vdw)
        if score > best[1]:
            best = (subset, score)
    return best[0] if best[0] is not None else candidates

# ===================== REVISED TS BUILDER =====================

# --- New helper: detect if adding a bond would create an undesired tiny cycle (e.g. triangle) ---
def _creates_small_cycle(base_graph: nx.Graph, i: int, j: int, max_cycle: int = 3) -> bool:
    if base_graph.has_edge(i, j):
        return False
    temp = base_graph.copy()
    temp.add_edge(i, j)
    try:
        cycles = nx.cycle_basis(temp)
    except Exception:
        return False
    for c in cycles:
        if i in c and j in c and len(c) <= max_cycle:
            return True
    return False

# --- New helper: filter candidates that would produce path shortcuts or tiny cycles ---
def _filter_vib_candidates(base_graph: nx.Graph,
                           candidates: List[Dict[str,Any]],
                           allow_shortcuts: bool = False,
                           avoid_small_cycles: bool = True,
                           small_cycle_size: int = 3) -> List[Dict[str,Any]]:
    filtered = []
    for c in candidates:
        i, j = c["bond"]
        reasons = []
        if c.get("already"):
            c["filter_reasons"] = ["existing_edge"]
            filtered.append(c)
            continue
        if not allow_shortcuts and nx.has_path(base_graph, i, j):
            try:
                sp_len = nx.shortest_path_length(base_graph, i, j)
                if sp_len <= 2:
                    reasons.append(f"shortest_path_{sp_len}")
            except Exception:
                pass
        if avoid_small_cycles and not reasons:
            if _creates_small_cycle(base_graph, i, j, max_cycle=small_cycle_size):
                reasons.append("tiny_cycle")
        if reasons:
            c["filter_reasons"] = reasons
            continue
        c["filter_reasons"] = []
        filtered.append(c)
    return filtered

# --- NEW: verbose TS candidate diagnostics helper ---
def _diagnose_vib_candidates(base_graph: nx.Graph,
                             frame_ts: Atoms,
                             raw_candidates: List[Dict[str,Any]],
                             filtered: List[Dict[str,Any]],
                             chosen_set: set,
                             enable_debug: bool):
    if not enable_debug:
        return
    filt_set = {tuple(c["bond"]) for c in filtered}
    for c in raw_candidates:
        i, j = c["bond"]
        path_len = None
        small_cycle_flag = False
        reason = []
        if nx.has_path(base_graph, i, j):
            try:
                path_len = nx.shortest_path_length(base_graph, i, j)
            except Exception:
                path_len = None
        # Recompute tiny cycle potential
        if not base_graph.has_edge(i, j):
            small_cycle_flag = _creates_small_cycle(base_graph, i, j, max_cycle=3)
        if tuple(c["bond"]) not in filt_set and not c.get("already"):
            if path_len is not None and path_len <= 2:
                reason.append(f"short_path={path_len}")
            if small_cycle_flag:
                reason.append("tiny_cycle")
        if c.get("already"):
            reason.append("pre_existing")
        if tuple(c["bond"]) in chosen_set:
            reason.append("SELECTED")
        logger.debug(f"TS vib cand {i}-{j}: ratio={c['ratio']:.3f} delta={c.get('delta',0):.3f} "
                     f"path={path_len} tiny_cycle={small_cycle_flag} "
                     f"{' '.join(reason) if reason else 'unfiltered'}")

# --- NEW constants for cyclic proton transfer detection ---
CPT_HETERO = {'O', 'N'}
CPT_DONOR_METALS = {'Mn','Fe','Co','Ni','Cu','Zn','Ru','Rh','Pd','Ir','Pt'}
CPT_MIN_DELTA_HETERO = 0.60   # Å minimal |Δ| for H–hetero to consider in cycle
CPT_MIN_DELTA_DONOR  = 0.45   # Å minimal |Δ| for donor–H leg
CPT_PRUNE_METAL_RATIO_MARGIN = 0.15  # metal–H ratio must be at least this smaller improvement vs hetero to keep
CPT_PRUNE_METAL_DELTA_FACTOR = 0.60  # prune metal–H if its delta < factor * hetero–H delta for same H
# --------------------------------------------------------------------------------

def _classify_cyclic_proton_transfer(frame_ts: Atoms,
                                     vib_bond_info: Dict[Tuple[int,int], Tuple[float,float]]
                                     ) -> Dict[str, Any]:
    """
    Identify cyclic proton transfer motif:
        Donor1–H_a ... H_a–Acceptor–H_b ... H_b–Donor2
    Returns dict with:
      forced_edges: set[(i,j)] edges to force include
      prune_edges:  set[(i,j)] metal–H edges deemed spurious
    """
    # Normalize bond key ordering
    def norm(a,b): return (a,b) if a < b else (b,a)
    # Build per-H partner list
    per_h = {}
    for (a,b),(delta,init) in vib_bond_info.items():
        if frame_ts[a].symbol == 'H' and frame_ts[b].symbol != 'H':
            per_h.setdefault(a, []).append((b, delta, init))
        elif frame_ts[b].symbol == 'H' and frame_ts[a].symbol != 'H':
            per_h.setdefault(b, []).append((a, delta, init))
    forced = set()
    to_prune = set()
    # For each hetero atom, collect attached hydrogens (with qualifying delta)
    # Also build reverse map: hetero -> hydrogens
    hetero_map = {}
    for h, partners in per_h.items():
        for (p, delta, init) in partners:
            if frame_ts[p].symbol in CPT_HETERO and abs(delta) >= CPT_MIN_DELTA_HETERO:
                hetero_map.setdefault(p, []).append((h, delta))
    # Detect cycles
    for acceptor, h_list in hetero_map.items():
        if len(h_list) < 2:
            continue
        # Try all H pair combinations
        for i in range(len(h_list)):
            for j in range(i+1, len(h_list)):
                h_a, d_a = h_list[i]
                h_b, d_b = h_list[j]
                # For each hydrogen find donor (metal or hetero donor like N) distinct from acceptor
                def find_donor(h):
                    best = None
                    for (partner, delta, init) in per_h.get(h, []):
                        if partner == acceptor:
                            continue
                        sym = frame_ts[partner].symbol
                        if (sym in CPT_DONOR_METALS or sym in CPT_HETERO) and abs(delta) >= CPT_MIN_DELTA_DONOR:
                            # prefer largest delta
                            if best is None or abs(delta) > abs(best[1]):
                                best = (partner, delta)
                    return best
                donor_a = find_donor(h_a)
                donor_b = find_donor(h_b)
                if donor_a and donor_b:
                    # We have a cyclic proton transfer motif
                    for edge in [(donor_a[0], h_a), (h_a, acceptor),
                                 (acceptor, h_b), (h_b, donor_b[0])]:
                        forced.add(norm(*edge))
                    # Optionally connect donors via acceptor path (not forcing direct donor–donor bond)
    # Prune spurious metal–H relative to stronger hetero–H for same H
    for h, partners in per_h.items():
        # Separate hetero vs metal
        hetero_partners = [(p,delta,init) for (p,delta,init) in partners if frame_ts[p].symbol in CPT_HETERO]
        if not hetero_partners:
            continue
        # Use strongest hetero
        strong_het = max(hetero_partners, key=lambda x: abs(x[1]))
        for (p, delta, init) in partners:
            sym = frame_ts[p].symbol
            if sym in CPT_DONOR_METALS:
                if abs(delta) < CPT_PRUNE_METAL_DELTA_FACTOR * abs(strong_het[1]):
                    to_prune.add(norm(h,p))
                else:
                    # Distance ratio comparison
                    dist_hp = frame_ts.get_distance(h,p)
                    dist_hh = frame_ts.get_distance(h,strong_het[0])
                    # if metal-H much longer than hetero-H by margin
                    if dist_hp - dist_hh > CPT_PRUNE_METAL_RATIO_MARGIN:
                        to_prune.add(norm(h,p))
    return {"forced_edges": forced, "prune_edges": to_prune}

def _apply_cyclic_proton_transfer_logic(frame_ts: Atoms,
                                        vib_bonds: List[Tuple[int,int]],
                                        vib_bond_info: Dict[Tuple[int,int], Tuple[float,float]]):
    """
    Modify vib_bonds list in-place:
      - remove pruned metal–H edges
      - add (or retain) forced cycle edges
    Returns (updated_vib_bonds, forced_set)
    """
    classification = _classify_cyclic_proton_transfer(frame_ts, vib_bond_info)
    forced = classification["forced_edges"]
    pruned = classification["prune_edges"]
    # Rebuild vib bond list without pruned
    vb_norm = set((min(a,b), max(a,b)) for a,b in vib_bonds)
    vb_norm -= pruned
    # Add forced
    vb_norm |= forced
    updated = sorted(vb_norm)
    return updated, forced

def build_ts_reference_graph(frame_ts: Atoms,
                             vib_bonds: List[Tuple[int,int]],
                             vib_bond_info: Dict[Tuple[int,int], Tuple[float, float]],
                             frames_displaced: List[Atoms],
                             method: str = 'cheminf',
                             charge: int = 0,
                             multiplicity: Optional[int] = None,
                             enable_scoring: bool = True,
                             allow_shortcuts: bool = False,
                             avoid_small_cycles: bool = True) -> nx.Graph:
    logger.debug(f"TS: building base graph (scoring={enable_scoring})")
    g = build_graph(frame_ts, method=method, charge=charge, multiplicity=multiplicity)
    for i in range(len(frame_ts)):
        if i in g.nodes and g.nodes[i].get('symbol') != frame_ts[i].symbol:
            g.nodes[i]['symbol'] = frame_ts[i].symbol
    if not vib_bonds:
        return g
    # --- NEW: cyclic proton transfer enforcement / pruning BEFORE candidate building ---
    vib_bonds, forced_cycle = _apply_cyclic_proton_transfer_logic(frame_ts, vib_bond_info=vib_bond_info, vib_bonds=vib_bonds)
    # Prepare standard candidates
    vdw = {frame_ts[i].symbol: _vdw_radius(frame_ts[i].symbol) for i in range(len(frame_ts))}
    raw_candidates = _prepare_vib_candidates(g, frame_ts, vib_bonds, vdw, vib_bond_info)
    # Mark forced-cycle candidates
    for c in raw_candidates:
        a,b = c["bond"]
        if (min(a,b), max(a,b)) in forced_cycle:
            c["forced_cycle"] = True
    filtered = _filter_vib_candidates(g, raw_candidates,
                                      allow_shortcuts=allow_shortcuts,
                                      avoid_small_cycles=avoid_small_cycles,
                                      small_cycle_size=3)
    if enable_scoring:
        chosen = _select_vib_subset(g, frame_ts, filtered, vdw, enable_scoring=True)
    else:
        chosen = filtered
    # ALWAYS include forced cycle edges
    chosen_set = {tuple(c["bond"]) for c in chosen}
    for c in filtered:
        if c.get("forced_cycle") and tuple(c["bond"]) not in chosen_set:
            chosen.append(c)
            chosen_set.add(tuple(c["bond"]))
    _diagnose_vib_candidates(g, frame_ts, raw_candidates, filtered, chosen_set, enable_scoring)
    added = 0
    annotated = 0
    for c in filtered:
        i, j = c["bond"]
        delta, _init = vib_bond_info.get((i,j), vib_bond_info.get((j,i),(0.0,0.0)))
        selected = (tuple(c["bond"]) in chosen_set)
        if not selected:
            continue
        if g.has_edge(i, j):
            g[i][j]['vib_identified'] = True
            g[i][j]['vib_delta'] = delta
            g[i][j]['vib_selected'] = True
            g[i][j]['vib_added_in_ts'] = False
            if c.get("forced_cycle"):
                g[i][j]['vib_forced_cycle'] = True
            g[i][j]['TS'] = True
            annotated += 1
            continue
        ratio = c.get("ratio", 1.0)
        dist = c.get("dist", frame_ts.get_distance(i, j))
        bo_candidates = [c.get("bo", 1.0)]
        if bo_candidates[0] < 1.0: bo_candidates.append(1.0)
        if ratio < 0.50: bo_candidates.append(1.5)
        success = False
        for bo in bo_candidates:
            try:
                g.add_edge(i, j,
                           bond_order=float(bo),
                           distance=dist,
                           bond_type=(frame_ts[i].symbol, frame_ts[j].symbol),
                           vib_identified=True,
                           vib_delta=delta,
                           vib_selected=True,
                           vib_added_in_ts=True,
                           TS=True,
                           vib_forced_cycle=bool(c.get("forced_cycle")))
                success = True
                added += 1
                break
            except Exception as e:
                logger.debug(f"TS vib add retry {i}-{j} bo={bo}: {e}")
        if not success:
            logger.debug(f"TS vib bond {i}-{j} skipped after retries")
    logger.debug(f"TS: vib candidates={len(raw_candidates)}, filtered={len(filtered)}, added={added}, annotated existing={annotated}, forced_cycle={len(forced_cycle)}")
    rings = find_rings(g)
    g.graph['base_ring_count'] = len(rings)
    g.graph['vib_new_rings'] = []
    return g

# ===================== REVISED DISPLACED BUILDER =====================

def build_displaced_graph(frame: Atoms,
                         ts_graph: nx.Graph,
                         vib_bonds: List[Tuple[int,int]],
                         vib_bond_info: Dict[Tuple[int,int], Tuple[float, float]],
                         frame_label: str,
                         method: str = 'cheminf',
                         charge: int = 0,
                         multiplicity: Optional[int] = None,
                         debug_ascii: bool = False,
                         use_ts_template: bool = False,
                         strategy: str = 'template',
                         removal_margin: float = 0.12,
                         keep_slack: float = 0.08,
                         recovery_rank: int = 2) -> nx.Graph:
    """
    Hysteretic pruning:
      keep if ratio <= thr + keep_slack
      remove if ratio > thr + removal_margin
      borderline (thr+keep_slack < ratio <= thr+removal_margin) retained & tagged vib_borderline=True

    After both displaced graphs are built (handled in analysis), recovery may re-add top-N delta vib bonds
    that were removed in BOTH frames (see recovery in analyze_displacement_graphs).
    """
    g = ts_graph.copy()
    if not vib_bonds:
        return g

    vdw = {frame[k].symbol: _vdw_radius(frame[k].symbol) for k in range(len(frame))}
    vib_edges_ts = {
        (min(i, j), max(i, j))
        for i, j, d in ts_graph.edges(data=True)
        if d.get('vib_identified', False)
    }

    removed = 0
    retained = 0
    borderline = 0

    for (i, j) in list(g.edges()):
        key = (min(i, j), max(i, j))
        if key not in vib_edges_ts:
            continue  # only prune TS vib bonds
        # Geometry in displaced frame
        if i >= len(frame) or j >= len(frame):
            continue
        ratio, dist, _ = _vdw_ratio(frame, i, j, vdw)
        thr = _bond_threshold(frame[i].symbol, frame[j].symbol)
        upper_keep = thr + keep_slack
        upper_remove = thr + removal_margin
        decision = None
        if ratio <= upper_keep:
            decision = "keep"
        elif ratio > upper_remove:
            g.remove_edge(i, j)
            removed += 1
            decision = "remove"
        else:
            # borderline keep
            if 'distance' in g[i][j]:
                g[i][j]['distance'] = dist
            else:
                g.add_edge(i, j, distance=dist)
            g[i][j]['vib_retained'] = True
            g[i][j]['vib_borderline'] = True
            g[i][j]['vib_retained_ratio'] = ratio
            g[i][j]['vib_threshold'] = thr
            g[i][j]['vib_prune_reason'] = f"borderline thr={thr:.3f} ratio={ratio:.3f}"
            borderline += 1
            retained += 1
            continue
        if decision == "keep":
            if 'distance' in g[i][j]:
                g[i][j]['distance'] = dist
            g[i][j]['vib_retained'] = True
            g[i][j]['vib_retained_ratio'] = ratio
            g[i][j]['vib_threshold'] = thr
            g[i][j]['vib_prune_reason'] = f"kept ratio={ratio:.3f}<=thr+{keep_slack:.2f}"
            retained += 1
        elif decision == "remove":
            logger.debug(f"{frame_label}: remove vib {i}-{j} ratio={ratio:.3f} thr={thr:.3f} (>{thr+removal_margin:.3f})")

    g.graph['vib_retained_count'] = retained
    g.graph['vib_removed_count'] = removed
    g.graph['vib_borderline_count'] = borderline
    logger.debug(f"{frame_label}: vib bonds retained={retained} (borderline={borderline}), removed={removed}")

    if debug_ascii and vib_edges_ts:
        try:
            focus_atoms = sorted({a for pair in vib_edges_ts for a in pair})
            expanded = _expand_with_neighbors(g, focus_atoms, depth=1)
            sub = g.subgraph(expanded).copy()
        except Exception as e:
            logger.warning(f"{frame_label}: ASCII debug failed: {e}")

    return g

# --- NEW: enforce removal of static vib bonds ---
def _enforce_static_vib_bond_prune(ts_graph: nx.Graph,
                                   g1: nx.Graph,
                                   g2: nx.Graph,
                                   variation_min: float = 0.06) -> Dict[str,int]:
    """
    Remove TS vib bonds from displaced graphs if the vibrational mode does NOT
    modulate them (ratio span below variation_min). Focuses on relative change,
    not absolute threshold classification.
    """
    rem1 = rem2 = 0
    for i, j, data in ts_graph.edges(data=True):
        if not data.get('vib_identified'):
            continue
        span = data.get('vib_ratio_span')
        if span is None:
            continue
        sym_i = ts_graph.nodes[i].get('symbol')
        sym_j = ts_graph.nodes[j].get('symbol')
        has_h = ('H' in (sym_i, sym_j))
        metal_set = {'Mn','Fe','Co','Ni','Cu','Zn','Ru','Rh','Pd','Ir','Pt'}
        hetero = {'O','N','S','P'}
        is_metal_h = has_h and ((sym_i in metal_set) or (sym_j in metal_set))
        involves_hetero = (sym_i in hetero) or (sym_j in hetero)
        base = variation_min
        if is_metal_h:
            span_thr = base * 1.25
        elif involves_hetero and has_h:
            span_thr = base * 0.85
        elif data.get('vib_forced_cycle'):
            span_thr = base * 0.80
        else:
            span_thr = base
        data['vib_static_span_threshold'] = span_thr  # record
        data['vib_static_span'] = span
        if span < span_thr:
            if g1.has_edge(i, j):
                g1.remove_edge(i, j); rem1 += 1
            if g2.has_edge(i, j):
                g2.remove_edge(i, j); rem2 += 1
            data['vib_static_removed'] = True
    return {"static_removed_frame1": rem1, "static_removed_frame2": rem2}

# --- NEW: recovery + detailed final report ---
def _recover_removed_vib_bonds(ts_graph: nx.Graph,
                               g1: nx.Graph,
                               g2: nx.Graph,
                               frame1: Atoms,
                               frame2: Atoms,
                               vib_bonds: List[Tuple[int,int]],
                               vib_bond_info: Dict[Tuple[int,int], Tuple[float,float]],
                               max_rank: int = 2):
    if not vib_bonds:
        return
    # Rank vib bonds by |delta|
    ranked = sorted(vib_bonds,
                    key=lambda b: abs(vib_bond_info.get(b, vib_bond_info.get((b[1], b[0]), (0,0)))[0]),
                    reverse=True)
    top = set(ranked[:max_rank])
    vdw = {sym: _vdw_radius(sym) for sym in {a.symbol for a in frame1} | {a.symbol for a in frame2}}
    def geom(frame, i, j):
        r, d, thr = _vdw_ratio(frame, i, j, vdw)
        thr = _bond_threshold(frame[i].symbol, frame[j].symbol)
        return r, d, thr
    recovered = []
    for (i, j) in top:
        in1 = g1.has_edge(i, j)
        in2 = g2.has_edge(i, j)
        # Only consider if vib_identified in TS and absent in BOTH
        if not ts_graph.has_edge(i, j) or not ts_graph[i][j].get('vib_identified', False):
            continue
        if in1 or in2:
            continue
        r1, d1, thr1 = geom(frame1, i, j)
        r2, d2, thr2 = geom(frame2, i, j)
        pick = None
        if r1 < r2:
            pick = (g1, r1, d1, thr1, "frame1")
        else:
            pick = (g2, r2, d2, thr2, "frame2")
        g_target, rr, dd, thr, label = pick
        # Add with neutral bond order 1.0
        g_target.add_edge(i, j,
                          bond_order=1.0,
                          distance=dd,
                          vib_recovered=True,
                          vib_retained=True,
                          vib_retained_ratio=rr,
                          vib_threshold=thr,
                          vib_prune_reason="recovered_absent_both",
                          TS=True)  # keep TS tag for recovered vib TS bond
        recovered.append(((i, j), label, rr, thr))
    if recovered:
        for (bond, label, rr, thr) in recovered:
            logger.info(f"Recovered vib bond {bond} into {label} (ratio={rr:.3f} thr={thr:.3f})")

def _debug_vib_summary(ts_graph: nx.Graph,
                       g1: nx.Graph,
                       g2: nx.Graph,
                       frame1: Atoms,
                       frame2: Atoms,
                       vib_bonds: List[Tuple[int,int]]):
    if not vib_bonds:
        return
    logger.info("VIB BOND GEOMETRY SUMMARY:")
    vdw = {sym: _vdw_radius(sym) for sym in {a.symbol for a in frame1} | {a.symbol for a in frame2}}
    for (i, j) in vib_bonds:
        line = f"{i}-{j}: "
        # TS presence
        ts_present = ts_graph.has_edge(i, j) and ts_graph[i][j].get('vib_identified', False)
        line += f"TS={'Y' if ts_present else 'N'} "
        for tag, frame, g in (('f1', frame1, g1), ('f2', frame2, g2)):
            r, d, _ = _vdw_ratio(frame, i, j, vdw)
            thr = _bond_threshold(frame[i].symbol, frame[j].symbol)
            present = g.has_edge(i, j)
            reason = g[i][j].get('vib_prune_reason') if present else 'removed'
            line += f"{tag}:(ratio={r:.3f} thr={thr:.2f} {'P' if present else 'X'} {reason}) "
        logger.debug(line)

# --- NEW: structured vib bond summary ---
def summarize_vib_bonds(ts_graph: nx.Graph) -> List[Dict[str,Any]]:
    """
    Structured per-bond summary for debugging / export.
    """
    rows = []
    for i, j, d in ts_graph.edges(data=True):
        if not d.get('vib_identified'):
            continue
        rows.append({
            "bond": (i, j),
            "symbols": (ts_graph.nodes[i].get('symbol'), ts_graph.nodes[j].get('symbol')),
            "delta": d.get('vib_delta'),
            "forced_cycle": d.get('vib_forced_cycle', False),
            "selected": d.get('vib_selected', False),
            "static_removed": d.get('vib_static_removed', False),
            "span": d.get('vib_ratio_span'),
            "span_thr": d.get('vib_static_span_threshold'),
            "state": d.get('vib_state'),
            "filter_reasons": d.get('vib_filter_reasons', []),
        })
    return rows

# ===================== ANALYSIS (patch only where displaced graphs built) =====================

# --- Adjust analyze_displacement_graphs: remove forced presence enforcement ---
def analyze_displacement_graphs(frames: List[Atoms],
                                frame_indices: List[int],
                                method: str = 'cheminf',
                                charge: int = 0,
                                multiplicity: Optional[int] = None,
                                internal_changes: Optional[Dict[str, Any]] = None,
                                ascii_neighbor_shells: int = 1,
                                ascii_scale: float = 3.0,
                                ascii_include_h: bool = True,
                                debug: bool = False
                                ) -> Dict[str, Any]:
    if len(frame_indices) < 2:
        logger.warning("Need at least 2 frames for displacement analysis")
        return {}
    ts_idx = internal_changes.get('ts_frame', 0) if internal_changes else 0
    f1_idx, f2_idx = frame_indices[:2]
    frame_ts = frames[ts_idx]; frame1 = frames[f1_idx]; frame2 = frames[f2_idx]
    vib_bonds = []; vib_bond_info = {}
    if internal_changes and 'bond_changes' in internal_changes:
        vib_bond_info = internal_changes['bond_changes']
        vib_bonds = list(vib_bond_info.keys())

    g_ts = build_ts_reference_graph(frame_ts, vib_bonds, vib_bond_info,
                                    frames_displaced=[frame1, frame2],
                                    method=method, charge=charge,
                                    multiplicity=multiplicity,
                                    enable_scoring=True,
                                    allow_shortcuts=False,
                                    avoid_small_cycles=True)

    _annotate_vib_bond_states(g_ts, frames, ts_idx, [f1_idx, f2_idx], vib_bonds)

    g1 = build_displaced_graph(frame1, g_ts, vib_bonds, vib_bond_info,
                               frame_label=f"frame{f1_idx}",
                               method=method, charge=charge,
                               multiplicity=multiplicity,
                               debug_ascii=debug,
                               strategy='template')
    g2 = build_displaced_graph(frame2, g_ts, vib_bonds, vib_bond_info,
                               frame_label=f"frame{f2_idx}",
                               method=method, charge=charge,
                               multiplicity=multiplicity,
                               debug_ascii=debug,
                               strategy='template')

    # NEW: enforce removal of static vib bonds
    enforce_stats = _enforce_static_vib_bond_prune(g_ts, g1, g2, variation_min=0.06)
    if debug and (enforce_stats["static_removed_frame1"] or enforce_stats["static_removed_frame2"]):
        logger.debug(f"Static vib bond pruning (low variation): {enforce_stats}")

    # NEW: recovery phase (optional)
    if debug:
        _recover_removed_vib_bonds(g_ts, g1, g2, frame1, frame2, vib_bonds, vib_bond_info, max_rank=2)
        _debug_vib_summary(g_ts, g1, g2, frame1, frame2, vib_bonds)

    comp = compare_graphs(g1, g2, frame1, frame2)
    comp['bonds_formed_symbols'] = [(frame1[i].symbol, frame1[j].symbol) for (i,j) in comp['bonds_formed']]
    comp['bonds_broken_symbols'] = [(frame1[i].symbol, frame1[j].symbol) for (i,j) in comp['bonds_broken']]
    vib_edges_in_ts = [(i,j,d) for i,j,d in g_ts.edges(data=True) if d.get('vib_identified', False)]
    vib_state_counts = g_ts.graph.get('vib_state_counts', {})

    results = {
        "ts_idx": ts_idx,
        "frame1_idx": f1_idx,
        "frame2_idx": f2_idx,
        "ts_stats": {
            "n_atoms": len(frame_ts),
            "n_bonds": g_ts.number_of_edges(),
            "n_rings_total": len(find_rings(g_ts)),
            "n_rings_base": g_ts.graph.get('base_ring_count', 0),
            "n_rings_from_vib": 0,
            "vib_new_rings": [],
            "n_fragments": nx.number_connected_components(g_ts),
            "vib_bonds_total": len(vib_bonds),
            "vib_bonds_added": len(vib_edges_in_ts),
            "vib_state_counts": vib_state_counts
        },
        "frame1_stats": {
            "n_atoms": len(frame1),
            "n_bonds": g1.number_of_edges(),
            "n_rings": len(find_rings(g1)),
            "n_fragments": nx.number_connected_components(g1),
            "vib_retained_bonds": g1.graph.get('vib_retained_count', 0),
            "vib_removed_bonds": g1.graph.get('vib_removed_count', 0),
            "vib_static_removed": enforce_stats["static_removed_frame1"],
        },
        "frame2_stats": {
            "n_atoms": len(frame2),
            "n_bonds": g2.number_of_edges(),
            "n_rings": len(find_rings(g2)),
            "n_fragments": nx.number_connected_components(g2),
            "vib_retained_bonds": g2.graph.get('vib_retained_count', 0),
            "vib_removed_bonds": g2.graph.get('vib_removed_count', 0),
            "vib_static_removed": enforce_stats["static_removed_frame2"],
        },
        "comparison": comp
    }
    # NEW: attach debug vib summaries
    if debug:
        try:
            results["ts_vib_overview"] = debug_ts_logic_overview(g_ts)
            results["vib_bond_table"] = summarize_vib_bonds(g_ts)
        except Exception as e:
            logger.debug(f"Vib overview attach failed: {e}")
    try:
        ascii_pkg = generate_transformation_ascii_with_ts(
            g_ts, g1, g2, comp,
            neighbor_shells=ascii_neighbor_shells,
            scale=ascii_scale,
            include_h=ascii_include_h
        )
        results.update(ascii_pkg)
    except Exception as e:
        logger.warning(f"ASCII generation failed: {e}", exc_info=debug)
    return results

def _annotate_vib_bond_states(ts_graph: nx.Graph,
                              frames: List[Atoms],
                              ts_idx: int,
                              displaced_indices: List[int],
                              vib_bonds: List[Tuple[int,int]]):
    """
    Populate per-bond vibrational presence metadata on TS graph only.
    Does not add/remove edges; purely annotates.
    """
    if not vib_bonds or len(displaced_indices) < 2:
        return
    vdw = {sym: _vdw_radius(sym) for sym in {a.symbol for a in frames[displaced_indices[0]]} | {a.symbol for a in frames[displaced_indices[1]]}}
    f1_idx, f2_idx = displaced_indices[:2]
    frame1 = frames[f1_idx]
    frame2 = frames[f2_idx]

    def ratio_present(frame: Atoms, i: int, j: int) -> Tuple[float,bool,float]:
        sym1, sym2 = frame[i].symbol, frame[j].symbol
        dist = frame.get_distance(i, j)
        r_sum = vdw.get(sym1, 2.0) + vdw.get(sym2, 2.0)
        ratio = dist / r_sum if r_sum > 1e-6 else 2.0
        thr = _bond_threshold(sym1, sym2)
        return ratio, (ratio < thr), thr

    state_counts = {"forming":0, "breaking":0, "persistent":0, "absent":0}

    for (a, b, data) in ts_graph.edges(data=True):
        if not data.get('vib_identified', False):
            continue
        i, j = a, b
        r1, p1, thr = ratio_present(frame1, i, j)
        r2, p2, _   = ratio_present(frame2, i, j)

        if p1 and not p2:
            vib_state = "breaking"
        elif not p1 and p2:
            vib_state = "forming"
        elif p1 and p2:
            vib_state = "persistent"
        else:
            vib_state = "absent"

        state_counts[vib_state] += 1
        span = abs(r1 - r2)  # NEW: variation metric across displaced frames

        ts_graph[i][j]['vib_frame_presence'] = {f1_idx: p1, f2_idx: p2}
        ts_graph[i][j]['vib_frame_ratio'] = {f1_idx: r1, f2_idx: r2}
        ts_graph[i][j]['vib_state'] = vib_state
        ts_graph[i][j]['vib_threshold'] = thr
        ts_graph[i][j]['vib_ratio_span'] = span  # NEW

    ts_graph.graph['vib_state_counts'] = state_counts
    logger.debug(f"TS vib bond states: {state_counts}")

def generate_transformation_ascii_with_ts(graph_ts: nx.Graph,
                                         graph_1: nx.Graph,
                                         graph_2: nx.Graph,
                                         comparison: Dict[str, Any],
                                         neighbor_shells: int = 1,
                                         scale: float = 3.0,
                                         include_h: bool = True) -> Dict[str, Any]:
    """
    Generate ASCII with TS as orientation reference.
    Includes atoms participating in:
      - Formed/broken/order-changed bonds (comparison dict)
      - All TS vibrational bonds (edges with vib_identified)
    Expanded by neighbor_shells in the TS graph.
    """
    # Core from formed/broken/order-changed
    core_from_comparison = set(_collect_transformation_nodes(comparison))
    # Core from TS vibrational bonds
    core_from_ts_vib = set()
    for i, j, d in graph_ts.edges(data=True):
        if d.get('vib_identified'):
            core_from_ts_vib.add(i); core_from_ts_vib.add(j)
    core = sorted(core_from_comparison | core_from_ts_vib)
    if not core:
        return {
            "ascii_ts": "<no change>",
            "ascii_ref": "<no change>",
            "ascii_disp": "<no change>",
            "nodes": [],
            "core_nodes": []
        }
    expanded = _expand_with_neighbors(graph_ts, core, depth=neighbor_shells)
    sub_ts = graph_ts.subgraph(expanded).copy()
    sub_1 = graph_1.subgraph(expanded).copy()
    sub_2 = graph_2.subgraph(expanded).copy()
    try:
        ascii_ts = graph_to_ascii(sub_ts, scale=scale, include_h=include_h, reference=None)
        ascii_1 = graph_to_ascii(sub_1, scale=scale, include_h=include_h, reference=sub_ts)
        ascii_2 = graph_to_ascii(sub_2, scale=scale, include_h=include_h, reference=sub_ts)
    except Exception as e:
        logger.warning(f"ASCII generation failed: {e}")
        return {
            "ascii_ts": "<error>",
            "ascii_ref": "<error>",
            "ascii_disp": "<error>",
            "nodes": expanded,
            "core_nodes": core
        }
    return {
        "ascii_ts": ascii_ts,
        "ascii_ref": ascii_1,
        "ascii_disp": ascii_2,
        "nodes": expanded,
        "core_nodes": core
    }


def compare_graphs(ref_graph: nx.Graph,
                   displaced_graph: nx.Graph,
                   frame_ref: Atoms,
                   frame_disp: Atoms) -> Dict[str, Any]:
    """Compare two molecular graphs and classify transformation."""
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
    
    # Fragment counts
    ref_fragments = nx.number_connected_components(ref_graph)
    disp_fragments = nx.number_connected_components(displaced_graph)
    
    transformation_type = classify_transformation(
        bonds_formed,
        bonds_broken,
        bond_order_changes,
        rings_formed,
        rings_broken,
        ref_fragments,
        disp_fragments
    )
    
    return {
        "bonds_formed": bonds_formed,
        "bonds_broken": bonds_broken,
        "bond_order_changes": bond_order_changes,
        "rings_formed": rings_formed,
        "rings_broken": rings_broken,
        "ref_fragments": ref_fragments,
        "disp_fragments": disp_fragments,
        "transformation_type": transformation_type,
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
    
    if len(bonds_formed) >= 2 and len(bonds_broken) >= 2:
        classifications.append("concerted_rearrangement")
    
    if (bonds_formed or bonds_broken) and ref_frags == disp_frags and not classifications:
        classifications.append("rearrangement")
    
    return " + ".join(classifications) if classifications else "no_significant_change"


def _collect_transformation_nodes(comparison: Dict[str, Any]) -> List[int]:
    """Collect unique atom indices involved in transformations"""
    involved = set()
    for i, j in comparison.get('bonds_formed', []):
        involved.update([i, j])
    for i, j in comparison.get('bonds_broken', []):
        involved.update([i, j])
    for (i, j) in comparison.get('bond_order_changes', {}).keys():
        involved.update([i, j])
    return sorted(involved)


def _expand_with_neighbors(graph: nx.Graph, core_nodes: List[int], depth: int = 1) -> List[int]:
    """Expand core nodes by neighbor shells"""
    expanded = set(core_nodes)
    frontier = set(core_nodes)
    for _ in range(depth):
        new_frontier = set()
        for n in frontier:
            new_frontier.update(graph.neighbors(n))
        new_frontier -= expanded
        expanded.update(new_frontier)
        frontier = new_frontier
        if not frontier:
            break
    return sorted(expanded)


def print_graph_analysis(graph_results: Dict[str, Any],
                        atom_index_map: Optional[Dict[int, str]] = None,
                        debug: bool = False):
    """Print formatted graph comparison results (summary; show raw ASCII blocks only)."""
    print("\n" + "="*80)
    print(" "*25 + "GRAPH ANALYSIS SUMMARY")
    print("="*80)

    if debug and 'ts_stats' in graph_results:
        ts_stats = graph_results['ts_stats']
        print(f"\nTS (frame {graph_results.get('ts_idx', 0)}):")
        print(f"  Atoms: {ts_stats['n_atoms']}, Bonds: {ts_stats['n_bonds']}")
        print(f"  Rings: {ts_stats.get('n_rings_total',0)} total ({ts_stats.get('n_rings_base',0)} base + {ts_stats.get('n_rings_from_vib',0)} from vib bonds)")
        print(f"  Fragments: {ts_stats['n_fragments']}")
        vib_total = ts_stats.get('vib_bonds_total', 0)
        vib_added = ts_stats.get('vib_bonds_added', 0)
        print(f"  Vib-active bonds: {vib_total} detected, {vib_added} added to TS")
        if vib_added < vib_total:
            print(f"    ({vib_total - vib_added} not added: already present or rejected)")

    if debug:
        f1_stats = graph_results['frame1_stats']
        f2_stats = graph_results['frame2_stats']
        print(f"\nDisplaced frame {graph_results['frame1_idx']}:")
        print(f"  Atoms: {f1_stats['n_atoms']}, Bonds: {f1_stats['n_bonds']}")
        print(f"  Rings: {f1_stats['n_rings']}, Fragments: {f1_stats['n_fragments']}")
        if any(f1_stats.get(k,0) for k in ('vib_retained_bonds','vib_removed_bonds')):
            print(f"  Vib bonds: retained={f1_stats.get('vib_retained_bonds',0)}, removed={f1_stats.get('vib_removed_bonds',0)}")
        print(f"\nDisplaced frame {graph_results['frame2_idx']}:")
        print(f"  Atoms: {f2_stats['n_atoms']}, Bonds: {f2_stats['n_bonds']}")
        print(f"  Rings: {f2_stats['n_rings']}, Fragments: {f2_stats['n_fragments']}")
        if any(f2_stats.get(k,0) for k in ('vib_retained_bonds','vib_removed_bonds')):
            print(f"  Vib bonds: retained={f2_stats.get('vib_retained_bonds',0)}, removed={f2_stats.get('vib_removed_bonds',0)}")

    comp = graph_results['comparison']
    print(f"\n{'─'*80}")
    print(f"Transformation: {comp['transformation_type']}")

    # ALWAYS show formed / broken bonds
    if comp['bonds_formed']:
        print(f"\nBonds Formed ({len(comp['bonds_formed'])}):")
        for (i, j), (sym1, sym2) in zip(comp['bonds_formed'], comp['bonds_formed_symbols']):
            print(f"  ({i}, {j}) [{sym1}-{sym2}]")
    if comp['bonds_broken']:
        print(f"\nBonds Broken ({len(comp['bonds_broken'])}):")
        for (i, j), (sym1, sym2) in zip(comp['bonds_broken'], comp['bonds_broken_symbols']):
            print(f"  ({i}, {j}) [{sym1}-{sym2}]")

    # Detailed sections remain debug-gated
    if debug and comp['bond_order_changes']:
        print(f"\nBond Order Changes ({len(comp['bond_order_changes'])}):")
        for (i, j), (old, new) in comp['bond_order_changes'].items():
            if atom_index_map:
                sym1, sym2 = atom_index_map[i], atom_index_map[j]
                print(f"  ({i}, {j}) [{sym1}-{sym2}]: {old:.1f} → {new:.1f}")
            else:
                print(f"  ({i}, {j}): {old:.1f} → {new:.1f}")
    if debug and comp.get('rings_formed'):
        print(f"\nRings Formed: {len(comp['rings_formed'])}")
    if debug and comp.get('rings_broken'):
        print(f"\nRings Broken: {len(comp['rings_broken'])}")

    if 'ascii_ts' in graph_results:
        print(f"\n{'ASCII TRANSFORMATION (TS → displaced)':=^80}")
        if debug:
            nodes = graph_results.get('nodes', [])
            core = graph_results.get('core_nodes', [])
            if nodes:
                print(f"Subgraph nodes ({len(nodes)}): {nodes}")
            if core:
                print(f"Core (changing) nodes ({len(core)}): {core}")
        print("\nTS (reference):\n")
        print(graph_results['ascii_ts'])
        print(f"\nDisplaced frame {graph_results['frame1_idx']}:\n")
        print(graph_results.get('ascii_ref', '<no ascii_ref>'))
        print(f"\nDisplaced frame {graph_results['frame2_idx']}:\n")
        print(graph_results.get('ascii_disp', '<no ascii_disp>'))
    print("="*80)

# === PIPELINE OVERVIEW (added) =================================================
# 1. TS base graph: build_graph() infers baseline bonding from TS Cartesian geometry.
# 2. Vibrational change candidates: internal_changes['bond_changes'] supplies (i,j)->(delta,initial_distance).
# 3. Pre-processing:
#    - _apply_cyclic_proton_transfer_logic:
#         * prunes weak metal–H when a stronger hetero–H dominates same H
#         * forces inclusion of cyclic proton-transfer legs (donor–H–acceptor–H–donor)
# 4. Candidate evaluation:
#    - geometric screen (_filter_vib_candidates) removes tiny cycles / shortcuts
#    - scoring (_score_subset / _select_vib_subset) keeps a consistent minimal set
#    - forced cycle edges always restored
# 5. TS augmentation: selected (and forced) edges annotated with vib_* flags + TS=True
# 6. Displaced graphs (two frames around TS):
#    - start as copies of TS graph
#    - hysteretic pruning removes vib edges too long (ratio > thr+margin) but retains borderline
# 7. Mode-based refinement:
#    - _annotate_vib_bond_states records per-frame ratios + presence flags
#    - _enforce_static_vib_bond_prune eliminates vib edges whose ratio span (|r1-r2|) is too small (not modulated by mode)
# 8. Comparison: compare_graphs() classifies formed/broken between displaced frames only
# 9. Optional recovery (debug) may re-add top-delta vib bonds removed in both frames
#
# Interpretation rule-of-thumb now:
#   - A bond is considered “made” in frame2 vs frame1 if present only in frame2
#   - “Broken” if present only in frame1
#   - Static vib TS edges (little modulation) are suppressed before comparison to reduce noise
#
# This matches the conceptual model you described. Remaining “smart” parts:
#   * selective forcing of proton transfer cycles
#   * pruning of spurious metal–H when overshadowed by hetero–H
#   * variation-based (span) static filtering rather than absolute ratio cutoff
#   * optional recovery to avoid over-pruning aggressive geometry screens.
#
# TODO ideas (not implemented):
#   - Replace subset enumeration with a greedy marginal gain if candidate set grows
#   - Weight scoring by normalized delta / (initial distance) for scale invariance
#   - Track alternative classifications (e.g., multi-step vs concerted) via temporal ordering of ratio changes.

def debug_ts_logic_overview(graph_ts: nx.Graph) -> Dict[str, Any]:
    """
    Lightweight snapshot of TS vib logic decisions for external auditing.
    Returns:
      counts: dict of vib edge tag frequencies
      forced_cycle_edges: list of edges flagged vib_forced_cycle
      static_removed_edges: list of vib edges pruned later (if annotate already ran)
      span_stats: min/avg/max ratio span for vib edges
    """
    vib_edges = [(i,j,d) for i,j,d in graph_ts.edges(data=True) if d.get('vib_identified')]
    forced = []
    static_removed = []
    spans = []
    tags = {}
    for i,j,d in vib_edges:
        if d.get('vib_forced_cycle'):
            forced.append((i,j))
        if d.get('vib_static_removed'):
            static_removed.append((i,j))
        if 'vib_ratio_span' in d:
            spans.append(d['vib_ratio_span'])
        for k,v in d.items():
            if k.startswith('vib_') and v is True:
                tags[k] = tags.get(k,0)+1
    if spans:
        span_stats = {
            "min": float(min(spans)),
            "avg": float(sum(spans)/len(spans)),
            "max": float(max(spans)),
            "n": len(spans)
        }
    else:
        span_stats = {"min":None,"avg":None,"max":None,"n":0}
    return {
        "counts": tags,
        "forced_cycle_edges": forced,
        "static_removed_edges": static_removed,
        "span_stats": span_stats,
        "total_vib_edges": len(vib_edges)
    }

__all__ = [
    "build_ts_reference_graph",
    "build_displaced_graph",
    "analyze_displacement_graphs",
    "compare_graphs",
    "print_graph_analysis",
    "generate_transformation_ascii_with_ts",
    "debug_ts_logic_overview",
]
