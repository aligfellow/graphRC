# vib_analysis

Vibrational mode structural interpretation + lightweight graph transformation detection.

[![PyPI Downloads](https://static.pepy.tech/badge/vib-analysis)](https://pepy.tech/projects/vib-analysis)

## Installation
This can be installed via `pypi` with: 
```bash
pip install vib_analysis
```
Or locally by:
```bash
git clone https://github.com/aligfellow/vib_analysis.git
cd vib_analysis
pip install .
```

## Core Pipeline

1. Trajectory extraction  
   - Parsed from XYZ or QM output (cclib first; ORCA + orca_pltvib fallback).
2. Internal coordinate screening  
   - Significant bonds / angles / dihedrals via thresholds.
3. TS graph construction (`graph_compare.build_ts_reference_graph`)  
   - Base graph from `xyzgraph.build_graph`.  
   - Vibrational bond candidates from internal bond_changes.  
   - Pre-processing:
     - Cyclic proton transfer logic:
       - Force inclusion of donor–H–acceptor–H–donor legs (`vib_forced_cycle=True`).
       - Prune weak metal–H when overshadowed by stronger hetero–H contacts.
     - Geometric filtering removes tiny artificial 3‑cycles and path shortcuts (distance redundancy).
   - Subset scoring prefers:
     - Ring closure (5/6) formation
     - Strong geometry (ratio below threshold)
     - Large |Δ| displacement
     - Penalizes redundant near-parallel edges.
4. TS augmentation  
   - Selected edges tagged: `vib_identified`, `vib_selected`, `vib_delta`, optional `vib_forced_cycle`.
5. Displaced frame graphs (two frames around TS)  
   - Start as TS copy.  
   - Hysteretic pruning:
     - Keep: ratio ≤ threshold + keep_slack
     - Remove: ratio > threshold + removal_margin
     - Borderline retained + annotated.
6. Mode modulation refinement  
   - Per vib TS edge store per-frame ratios + `vib_ratio_span = |r1 - r2|`.
   - Static pruning: remove edges with insufficient span (< adaptive threshold):
     - Metal–H stricter (×1.25)
     - Hetero–H slightly looser (×0.85)
     - Forced cycle looser (×0.80).
   - Marks `vib_static_removed`.
7. Optional recovery (debug)  
   - Re-add top-|Δ| vib edges missing in both displaced graphs.
8. Comparison  
   - Formed / broken edges detected only between displaced graphs.
   - ASCII tri-panel (TS / frame1 / frame2) with consistent orientation.

## Key Edge Attributes

| Attribute | Meaning |
|-----------|---------|
| vib_identified | Came from vibrational candidate list |
| vib_selected | Survived scoring/filtering |
| vib_forced_cycle | Added by cyclic proton transfer inference |
| vib_delta | Displacement magnitude driver |
| vib_ratio_span | Variation between displaced frames |
| vib_static_removed | Pruned for low modulation |
| vib_retained / vib_borderline | Survived hysteretic pruning |

## CLI (excerpt)

```
python -m vib_analysis.cli input.out --mode 0 --graph --debug \
  --ascii-shells 1 --ascii-scale 2.5 --show-h
```

Important flags:
- `--graph` enable graph transformation layer
- `--frames i j` override auto-selected displaced frames
- `--ascii-shells N` neighbor expansion around transformation core (default 1)
- `--debug` attaches vib overview tables: `ts_vib_overview`, `vib_bond_table`
- `--displacement-amplitude k` saves symmetric displaced pair (k=1..4)

## Interpretation Guide

- Bond “formed” = present only in the second displaced graph.
- Bond “broken” = present only in the first displaced graph.
- Proton transfer cycles: look for contiguous donor–H–acceptor–H–donor path with all edges TS‑tagged.
- Static noise removed: edges with tiny span (< ~0.06 adjusted) suppressed before comparison.

## Debug Structures

When `--debug`:
- `ts_vib_overview`: counts of vib_* tags + span statistics.
- `vib_bond_table`: per-edge rows (bond, symbols, delta, span, thresholds, state).

## Dependencies

Core:
- xyzgraph
- ase
- networkx
Optional (input parsing):
- cclib
ORCA helper:
- orca + orca_pltvib in PATH

## Return Objects (API)

`run_vib_analysis(...)` returns:
```
{
  frames, freqs,
  vibrational: { bond_changes, angle_changes, ... },
  graph: { comparison, ts_stats, frame1_stats, frame2_stats, ascii_ts, ... },
  ts_vib_overview (debug),
  vib_bond_table (debug)
}
```

## Usage
>[!Note]
> Requires a `*trj.xyz` file of the structure, or output files from various QM softwares (see below)
```
[n_atoms]
comment line
<atomic symbol/number> <x> <y> <z>
... ... ... ... 
```

- **Auto-detection**: If input is `.xyz`, trajectory is read directly. Otherwise, attempts parsing with cclib first, then falls back to orca_pltvib.
- **Mode default**: `--mode` defaults to 0 (imaginary mode).
- Trajectory files are saved to disk when possible; if write fails, analysis proceeds without saving.

`--orca_path` can be provided, if absent this will default to checking for ORCA in PATH with `os.system("which orca")`

>[!IMPORTANT]
>- **atom indices are zero indexed** (though the viewer used below is *one indexed*)

### Future addition?
>[!TIP]  
> Suggestions?

## Command line interface
```bash
> vib_analysis -h
usage: vib_analysis [-h] [--mode MODE] [--orca_path ORCA_PATH] [--save-displacement]
                    [--no-save] [--bond_tolerance BOND_TOLERANCE] [--angle_tolerance ANGLE_TOLERANCE]
                    [--dihedral_tolerance DIHEDRAL_TOLERANCE] [--bond_threshold BOND_THRESHOLD]
                    [--angle_threshold ANGLE_THRESHOLD] [--dihedral_threshold DIHEDRAL_THRESHOLD]
                    [--ts_frame TS_FRAME] [--all] [--graph]
                    input

positional arguments:
  input                 Input file (XYZ trajectory, ORCA output, or Gaussian log)

options:
  -h, --help            show this help message and exit
  --mode MODE           Mode index to analyze (default: 0, zero-indexed)
  --orca_path ORCA_PATH
                        Path to ORCA binary (optional)
  --save-displacement, -sd
                        Save displaced structures (frames 1 and -1 from trajectory)
  --no-save             Do not save trajectory file to disk (keep in memory only)
  --bond_tolerance BOND_TOLERANCE
                        Bond detection tolerance multiplier. Default: 1.4
  --angle_tolerance ANGLE_TOLERANCE
                        Angle detection tolerance multiplier. Default: 1.1
  --dihedral_tolerance DIHEDRAL_TOLERANCE
                        Dihedral detection tolerance multiplier. Default: 1.0
  --bond_threshold BOND_THRESHOLD
                        Minimum bond change to report (Å). Default: 0.4
  --angle_threshold ANGLE_THRESHOLD
                        Minimum angle change to report (degrees). Default: 10
  --dihedral_threshold DIHEDRAL_THRESHOLD
                        Minimum dihedral change to report (degrees). Default: 20
  --ts_frame TS_FRAME   Reference frame index. Default: 0
  --all                 Report all internal coordinate changes
  --graph, -g           Perform graph analysis of structural transformations
```

## Atom Symbols in Output
Each reported internal coordinate includes element symbols:
```
Bond (11, 12) [C-O]: Δ = 1.432 Å, Initial Length = 2.064 Å
Angle (13, 12, 29) [C-N-H]: Δ = 11.020°, Initial Value = 122.116°
Dihedral (31, 13, 14, 32) [C-C-C-C]: Δ = 29.557°, Initial Value = 185.910°
```
The results dictionary contains:
```python
results['atom_index_map']  # { index: symbol }
```

## Python interface
See examplese/examples.ipynb
This function will return a dictionary of the results, and printing can be turned on to produce the same as the CLI
For example:
```python
from vib_analysis import run_vib_analysis

orca_out = 'data/bimp.v000.xyz'

results = run_vib_analysis(
        input_file=orca_out,
    )

print(results)

theoretical_bond_changes = [(11,12), (10,14)]
if all(bond in results['bond_changes'] for bond in theoretical_bond_changes):
    print(f'True: All theoretical bond changes {theoretical_bond_changes} found in results.')
```
Outputs:
```python
{'bond_changes': {(11, 12): (2.052, 2.064)}, 'angle_changes': {}, 'minor_angle_changes': {(13, 12, 29): (14.436, 122.116), (29, 12, 30): (12.54, 117.79)}, 'dihedral_changes': {(32, 14, 15, 20): (43.451, 350.826)}, 'minor_dihedral_changes': {(2, 1, 10, 11): (49.302, 194.336), (29, 12, 13, 31): (67.358, 17.521)}, 'frame_indices': [5, 15], 'atom_index_map': {0: 'O', 1: 'C', ... }}
True: All theoretical bond changes [(11, 12), (10, 14)] found in results.
```
  - This can be used to check for a known vibrational mode (theoretical_bond_change) in `results['bond_changes']`
  - So in theory this could identify whether the correct TS mode has been identidied in a high throughput search if the atom indices are known (or available automatically)

## More detailed information
- the `--all` flag turns on reporting of coupled internal coordinate changes, including:
   - Default output:
  ```bash
  =========================== Significant Bond Changes ===========================
  ========================== Significant Angle Changes ===========================
  ========================= Significant Dihedral Changes =========================
  ```
    - additional output - not necessarily insignificant changes in internal coordinates but strongly coupled
  ```bash
  ============================= Minor Angle Changes ==============================
  ============================ Minor Dihedral Changes ============================  
  ```
  - *i.e.* where a bond is changed, the angles around it will be altered across a vibrational trajectory and those angles would be significant enough to report as a change
  - where one of these atoms is involved in a *significant* bond change, the angle is classed as minor due to the coupled nature of the internal coordinates
     - same applies for dihedrals
   
### Save Displaced Structures
Export frames at small displacements:
```bash
vib_analysis input.out --save-displacement
# or
vib_analysis input.out -sd
# Creates: input_F.xyz, input_R.xyz (Forward/Reverse along mode)
```
Works even when trajectory is kept in-memory only (with `--no-save`).
- this is convenient for running a tight optimisations in a pseudo IRC, "quick" reaction coordinate

## Graph Analysis

Analyze structural transformations by comparing molecular graphs:

```bash
vib_analysis input.out --graph
# or
vib_analysis input.out -g
```

This feature identifies:
- **Bond formation/cleavage**: New bonds formed or broken
- **Bond order changes**: Single ↔ double ↔ triple transitions
- **Cyclization/ring-opening**: Ring formation or breaking
- **Pericyclic reactions**: Concerted bond reorganization
- **Rearrangements**: Structural changes without fragmentation
- **Fragmentations**: Molecule breaking into pieces

Example output:
```
================================================================================
                         GRAPH ANALYSIS SUMMARY
================================================================================

Reference Structure:
  Atoms: 42
  Bonds: 44
  Rings: 3
  Fragments: 1

────────────────────────────────────────────────────────────────────────────────
Frame 5 vs Reference:
  Transformation Type: bond_formation + cyclization

  Bonds Formed (2):
    (10, 14) [C-C]
    (11, 12) [O-C]

  Bond Order Changes (1):
    (13, 14) [C-C]: 1.0 → 2.0

  Rings Formed: 1
    Size 5: (10, 11, 12, 13, 14)
```

## Minimal Examples 
### Example 1
Sample python use in examples/ folder:
![sn2 imaginary mode](images/sn2.gif)
    - visualisation using [v.2.0](https://github.com/briling/v) by [**Ksenia Briling @briling**](https://github.com/briling) 
    - `v sn2.v000.xyz` press `f` and then `q` ; then ```bash convert -delay 5 -loop 0 sn2*xpm sn2.gif```

From the command line:
```bash
> vib_analysis sn2.v000.xyz
 # OR
> vib_analysis sn2.out --parse_orca --mode 0

Analysed vibrational trajectory from examples/data/sn2.v000.xyz:

=========================== Significant Bond Changes ===========================
Bond (0, 4)  [C-F]   Δ =   1.584 Å,  Initial =   1.717 Å
Bond (0, 5)  [C-Cl]  Δ =   1.355 Å,  Initial =   1.952 Å
```
The magnitude and change (Δ) of the modes is somewhat meaningless, though the initial value in the reference frame is also reported.

### Example 2
![dihedral imaginary mode](images/dihedral.gif)
```bash
> vib_analysis dihedral.v000.xyz
# OR
> vib_analysis dihedral.out # defaults parsing via CCLIB and orca_pltvib (if available) 

Analysed vibrational trajectory from examples/data/dihedral.v000.xyz:

========================= Significant Dihedral Changes =========================
Dihedral (6, 0, 3, 7)  [F-C-C-F]  Δ =  39.557 °,  Initial = 359.998 °
```

>[!NOTE]
>The bond changes are hierarchical, so an angle with a large change as a consequence of a bonding change is not reported as a *significant* change.

### Example 3
![larger molecule sn2](images/sn2_large.gif)
```bash
> vib_analysis sn2_large.v000.xyz

Analysed vibrational trajectory from examples/data/sn2_large.v000.xyz:

=========================== Significant Bond Changes ===========================
Bond (0, 21)  [C-N]  Δ =   2.388 Å,  Initial =   2.158 Å
Bond (0, 1)   [C-I]  Δ =   1.878 Å,  Initial =   2.563 Å
```

### Example 4 - more involved 

Complex transformation with BIMP catalysed rearrangement
- including the `--all` flag to print *all* internal coordinate changes
![bimp rearrangement](images/bimp.gif)
```bash
> vib_analysis bimp.v000.xyz --all 

Analysed vibrational trajectory from examples/data/bimp.v000.xyz:

=========================== Significant Bond Changes ===========================
Bond (11, 12)  [O-C]  Δ =   2.052 Å,  Initial =   2.064 Å
Bond (10, 14)  [C-C]  Δ =   0.426 Å,  Initial =   2.656 Å

============================= Minor Angle Changes ==============================
Angle (13, 12, 29)  [C-C-H]  Δ =  14.436 °,  Initial = 122.116 °
Angle (12, 13, 14)  [C-C-C]  Δ =  14.118 °,  Initial = 123.702 °
Angle (29, 12, 30)  [H-C-H]  Δ =  12.540 °,  Initial = 117.790 °

Note: These angles are dependent on other changes and may not be significant on their own.

============================ Minor Dihedral Changes ============================
Dihedral (29, 12, 13, 31)  [H-C-C-N]  Δ =  67.358 °,  Initial =  17.521 °
Dihedral (12, 13, 31, 33)  [C-C-N-S]  Δ =  62.151 °,  Initial = 330.369 °
Dihedral (4, 9, 10, 11)    [C-C-C-O]  Δ =  50.966 °,  Initial = 169.776 °
Dihedral (0, 1, 10, 11)    [O-C-C-O]  Δ =  36.480 °,  Initial =  14.986 °

Note: These dihedrals are dependent on other changes and may not be significant on their own.
```
![bimp rearrangement zoom](images/bimp_zoom.gif)

- correctly identifies the bond change between atoms 11 and 12
- also identifies the lower magnitude "looser" bonding change of 10 and 14

```bash
> vib_analysis examples/data/bimp.v000.xyz 

Analysed vibrational trajectory from examples/data/bimp.v000.xyz:

=========================== Significant Bond Changes ===========================
Bond (11, 12)  [O-C]  Δ =   2.052 Å,  Initial =   2.064 Å
Bond (10, 14)  [C-C]  Δ =   0.426 Å,  Initial =   2.656 Å
```

### Example 5 
Mn catalyst hydrogenation
![Mn hydrogenation](images/mn.gif)
```bash
> vib_analysis mn.log --all
Parsing mn.log with cclib...
Written trajectory to: mn.v000.xyz

First 5 non-zero vibrational frequencies:
  Mode 0: -748.48 cm**-1  (imaginary)
  Mode 1: 20.26 cm**-1 
  Mode 2: 25.12 cm**-1 
  Mode 3: 32.45 cm**-1 
  Mode 4: 36.68 cm**-1 

Analysed vibrational trajectory (Mode 0 with frequency -748.48 cm**-1):

=========================== Significant Bond Changes ===========================
Bond (5, 65)   [N-H]   Δ =   1.776 Å,  Initial =   1.319 Å
Bond (65, 66)  [H-O]   Δ =   1.665 Å,  Initial =   1.203 Å
Bond (64, 66)  [H-O]   Δ =   0.920 Å,  Initial =   1.711 Å
Bond (1, 65)   [Mn-H]  Δ =   0.875 Å,  Initial =   2.591 Å
Bond (1, 64)   [Mn-H]  Δ =   0.649 Å,  Initial =   1.898 Å

============================= Minor Angle Changes ==============================
Angle (64, 63, 1)  [H-H-Mn]  Δ =  45.388 °,  Initial =  88.897 °
Angle (5, 1, 63)   [N-Mn-H]  Δ =  16.471 °,  Initial =  96.799 °
Angle (61, 1, 63)  [C-Mn-H]  Δ =  15.528 °,  Initial =  81.202 °
Angle (2, 1, 63)   [P-Mn-H]  Δ =  13.032 °,  Initial = 171.266 °

Note: These angles are dependent on other changes and may not be significant on their own.
```
- this correctly identifies bonding changes of this transition state
- parsing the output prints the imaginary modes from the output file
- gaussian parsing with [cclib](https://github.com/cclib/cclib)

### ORCA parsing
Orca output parsing is also possible with `--parse_cclib` and separately with `--parse_orca` 
  - it appears that cclib cannot yet deal with orca_6.1.0 
```bash
> vib_analysis dihedral.out

Parsing dihedral.out with orca_pltvib...
INFO: Multiple 'VIBRATIONAL FREQUENCIES' sections found. Using the last one. # orca_pltvib by default runs on the FIRST occurance of a frequency section
Written trajectory to: vib_analysis/examples/data/dihedral.v000.xyz

First 5 non-zero vibrational frequencies:
  Mode 0: -225.09 cm**-1  (imaginary)
  Mode 1: 270.89 cm**-1 
  Mode 2: 612.80 cm**-1 
  Mode 3: 846.45 cm**-1 
  Mode 4: 899.93 cm**-1 

Analysed vibrational trajectory (Mode 0 with frequency -225.09 cm**-1):

========================= Significant Dihedral Changes =========================
Dihedral (6, 0, 3, 7)  [F-C-C-F]  Δ =  43.778 °,  Initial = 359.998 °
```

And again, with the bimp example:
```bash
vib_analysis bimp.out

Parsing dihedral.out with orca_pltvib...
INFO: Multiple 'VIBRATIONAL FREQUENCIES' sections found. Using the last one.
Written trajectory to: vib_analysis/examples/data/dihedral.v000.xyz

First 5 non-zero vibrational frequencies:
  Mode 0: -225.09 cm**-1  (imaginary)
  Mode 1: 270.89 cm**-1 
  Mode 2: 612.80 cm**-1 
  Mode 3: 846.45 cm**-1 
  Mode 4: 899.93 cm**-1 

Analysed vibrational trajectory (Mode 0 with frequency -225.09 cm**-1):

========================= Significant Dihedral Changes =========================
Dihedral (6, 0, 3, 7)  [F-C-C-F]  Δ =  43.778 °,  Initial = 359.998 °
```
- this output used `orca_6.0.1` and the `CCLIB` parsing is supported 
- newer versions, *i.e.* `orca_6.1.0`, will fall back to `orca_pltvib` if available
