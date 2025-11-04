# vib_analysis

> **Automated vibrational mode analysis and structural transformation detection**

Identify bond formation/breaking, angle changes, and dihedral rotations from vibrational trajectories with optional graph-based transformation analysis.

[![PyPI Downloads](https://static.pepy.tech/badge/vib-analysis)](https://pepy.tech/projects/vib-analysis)

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Examples](#examples)
- [Command Line Interface](#command-line-interface)
- [Python API](#python-api)
- [Advanced Options](#advanced-options)
- [Important Notes](#important-notes)

---

## Features

### Core Analysis
✅ **Automatic trajectory extraction** from XYZ files or QM output (ORCA, Gaussian via cclib)  
✅ **Internal coordinate tracking** - identifies significant bond, angle, and dihedral changes  
✅ **Smart filtering** - separates primary changes from coupled secondary effects  

### Advanced Analysis (--graph flag)
🔍 **Mode characterization** - identifies rotations, inversions, aromatic systems  
🔍 **Bond formation/cleavage detection**  
🔍 **Bond order changes** (single ↔ double ↔ triple)  
🔍 **Formal charge redistribution** tracking  
🔍 **ASCII molecular visualization** of transformations

>[!IMPORTANT]
> Bond orders and formal charges are **empirically assigned** by [xyzgraph](https://github.com/aligfellow/xyzgraph) and should be treated as **indicative only**.   
> They are particularly unreliable for metal-containing systems. Use them as qualitative guides, not quantitative predictions.  

---

## Installation

### From pypi
```bash
pip install vib-analysis
```

### From Source (*up-to-date*)
```bash
git clone https://github.com/aligfellow/vib_analysis.git
cd vib_analysis
pip install .
# or simply
pip install git+https://github.com/aligfellow/vib_analysis.git
```

### Dependencies
**Required:**
- `numpy` - Numerical operations
- `networkx` - Graph operations
- `xyzgraph` - Molecular graph construction (does the graph analysis)
- `cclib` - Parsing Gaussian/ORCA output

**Optional**:
- ORCA with `orca_pltvib` in PATH

---

## Quick Start

```bash
# Simple bond analysis
vib_analysis trajectory.xyz

# With analysis (characterization + graph + ASCII visualization)
vib_analysis calculation.out --graph

# Save structures for IRC calculations
vib_analysis calculation.out --save-displacement
```

---

## How It Works

### Key Components

**Core Analysis:**
- Selects relevant frames for comparison
- Identifies which bonds/angles/dihedrals change
- Detects coupled proton transfers with reduced threshold for H movements
- Compares graphs to detect transformations
- Filters and classifies changes

**xyzgraph's Role:**
- Constructs molecular graphs from 3D coordinates
- Assigns bond orders using empirical rules
- Calculates formal charges using valence rules
- Provides the graph infrastructure that we use

---

## Examples

> **Note:** All atom indices are **zero-indexed**

### Example 1: SN2 Reaction

![SN2 Animation](images/sn2.gif)
    - visualisation using [v.2.0](https://github.com/briling/v) by [**Ksenia Briling @briling**](https://github.com/briling):  
    - `v sn2.v000.xyz` press `f` and then `q` ; then ```bash convert -delay 5 -loop 0 sn2*xpm sn2.gif```

```bash
vib_analysis examples/data/sn2.v000.xyz
```

**Output:**
```
================================================================================
                              VIB_ANALYSIS
================================================================================

================================================================================
                    VIBRATIONAL TRAJECTORY ANALYSIS
================================================================================

=========================== Significant Bond Changes ===========================
Bond (0, 4)  [C-F]   Δ =   1.584 Å,  Initial =   1.717 Å
Bond (0, 5)  [C-Cl]  Δ =   1.355 Å,  Initial =   1.952 Å

================================================================================
```

**Interpretation:** Classic SN2 mechanism - concerted C-F bond breaking and C-Cl bond forming.

---

### Example 2: Dihedral Rotation

![Dihedral Rotation](images/dihedral.gif)

```bash
vib_analysis examples/data/dihedral.v000.xyz
```

**Output (truncated):**
```
================================================================================
                              VIB_ANALYSIS
================================================================================

================================================================================
                         MODE CHARACTERIZATION
================================================================================

Mode Type: ROTATION
Description: Single bond rotation

1 dihedral rotation(s) detected:
  (6, 0, 3, 7): Single bond C-C rotation (43.8°)

================================================================================
                    VIBRATIONAL TRAJECTORY ANALYSIS
================================================================================

========================= Significant Dihedral Changes =========================
Dihedral (6, 0, 3, 7)  [F-C-C-F]  Δ =  43.778 °,  Initial = 359.998 °

================================================================================
```

**Interpretation:** Internal rotation about C-C bond causing F-C-C-F dihedral change of ~44°.

---

### Example 3: Complex Rearrangement (Basic Analysis)

![BIMP Rearrangement](images/bimp.gif)

```bash
vib_analysis examples/data/bimp.v000.xyz
```

**Output:**
```
================================================================================
                              VIB_ANALYSIS
================================================================================

================================================================================
                    VIBRATIONAL TRAJECTORY ANALYSIS
================================================================================

=========================== Significant Bond Changes ===========================
Bond (11, 12)  [O-C]  Δ =   2.052 Å,  Initial =   2.064 Å
Bond (10, 14)  [C-C]  Δ =   0.426 Å,  Initial =   2.656 Å

================================================================================
```

**Interpretation:** Two significant bond changes detected - O-C formation and C-C breaking (formal \[2,3]-rearrangement).

---

### Example 4: With Graph Analysis & Charge Redistribution

![BIMP Rearrangement zoom](images/bimp_zoom.gif)

```bash
vib_analysis examples/data/bimp.out -g
```

**Output (excerpt):**
```
vib_analysis examples/data/bimp.out -g -as 2
================================================================================
                              VIB_ANALYSIS
================================================================================

Analyzed Mode 0: -333.88 cm⁻¹ (imaginary)

First 5 non-zero vibrational frequencies:
  Mode 0: -333.88 cm⁻¹ (imaginary)
  Mode 1: 8.57 cm⁻¹
  Mode 2: 12.72 cm⁻¹
  Mode 3: 13.27 cm⁻¹
  Mode 4: 15.83 cm⁻¹

================================================================================
                         VIBRATIONAL GRAPH ANALYSIS SUMMARY
================================================================================

Bonds Formed (1):
  Bond (11, 12) [O-C]: formed as single (order=1.0)

Bonds Broken (1):
  Bond (10, 14) [C-C]: broken from single (order=1.0)

Interpret with care, bond orders and charges are empirical and LOW confidence.

Bond Order Changes (4 bonds):
  Bond (1, 10) [C-C]: single→double (order 1.0→2.0)
  Bond (13, 14) [C-C]: single→double (order 1.0→2.0)
  Bond (1, 2) [C-N]: double→single (order 2.0→1.0)
  Bond (13, 31) [C-N]: double→single (order 2.0→1.0)

Formal Charge Redistribution (4 atoms):
  Atom 2 [N]: charge +1→+0 (Δq = -1.00)
  Atom 11 [O]: charge -1→+0 (Δq = +1.00)
  Atom 12 [C]: charge -1→+0 (Δq = +1.00)
  Atom 31 [N]: charge +1→+0 (Δq = -1.00)

================================================================================
ASCII REPRESENTATIONS
================================================================================

Transition State (TS):

           C‖
           ‖‖
           ‖‖
            ‖‖
            ‖‖    ------O
C-----------C‖----       **
            *              *
            *               **
           *                 /C
           *                /
          -C==           ///
      ---- =========    /
  ----        ========C/
C-                  ===

Frame 1:

           C
           |
           |
            |
            |     ------O
C-----------C-----
            |
            |
           |                 /C
           |                /
          -C--           ///
      ----    ------    /
  ----              --C/
C-

Frame 2:

           C‖
           ‖‖
           ‖‖
            ‖‖
            ‖‖    ------O
C-----------C‖----       \\
                           \
                            \\
                             /C
                            /
          -C==           ///
      ---- =========    /
  ----        ========C/
C-                  ===

================================================================================
                         MODE CHARACTERIZATION
================================================================================

Mode Type: BOND_CHANGE
Description: Bond formation/breaking

================================================================================
                    VIBRATIONAL TRAJECTORY ANALYSIS
================================================================================

=========================== Significant Bond Changes ===========================
Bond (11, 12)  [O-C]  Δ =   2.052 Å,  Initial =   2.064 Å
Bond (10, 14)  [C-C]  Δ =   0.426 Å,  Initial =   2.656 Å

================================================================================
```

**Interpretation:** Graph analysis reveals a rearrangement with bond formation/breaking, bond order changes, and charge redistribution. 

---

### Example 5: Showing All Changes (Including Minor)

```bash
vib_analysis examples/data/bimp.v000.xyz --all
```

Shows additional "Minor Angle Changes" and "Minor Dihedral Changes" sections with coupled secondary effects.

---

### Example 6: Larger SN2 System

![Large SN2](images/sn2_large.gif)

```bash
vib_analysis examples/data/sn2_large.v000.xyz
```

**Output:**
```
================================================================================
                              VIB_ANALYSIS
================================================================================

================================================================================
                    VIBRATIONAL TRAJECTORY ANALYSIS
================================================================================

=========================== Significant Bond Changes ===========================
Bond (0, 21)  [C-N]  Δ =   2.388 Å,  Initial =   2.158 Å
Bond (0, 1)   [C-I]  Δ =   1.878 Å,  Initial =   2.563 Å

================================================================================
```

**Interpretation:** SN2 reaction in larger molecular context - C-I bond breaking and C-N bond forming.

---

### Example 7: Mn Catalyst Hydrogenation

![Mn Hydrogenation](images/mn-h2.gif)

```bash
vib_analysis examples/data/mn-h2.log --all
```

**Output:**
```
================================================================================
                              VIB_ANALYSIS
================================================================================

Parsing mn-h2.log with cclib...
Saved trajectory to mn-h2.v000.xyz
Loaded 20 frames from trajectory

Analyzed Mode 0: -748.48 cm⁻¹ (imaginary)

First 5 non-zero vibrational frequencies:
  Mode 0: -748.48 cm⁻¹ (imaginary)
  Mode 1: 20.26 cm⁻¹
  Mode 2: 25.12 cm⁻¹
  Mode 3: 32.45 cm⁻¹
  Mode 4: 36.68 cm⁻¹

================================================================================
                    VIBRATIONAL TRAJECTORY ANALYSIS
================================================================================

=========================== Significant Bond Changes ===========================
Bond (5, 65)   [N-H]   Δ =   1.776 Å,  Initial =   1.319 Å
Bond (65, 66)  [H-O]   Δ =   1.665 Å,  Initial =   1.203 Å
Bond (64, 66)  [H-O]   Δ =   0.920 Å,  Initial =   1.711 Å
Bond (1, 64)   [Mn-H]  Δ =   0.649 Å,  Initial =   1.898 Å
Bond (63, 64)  [H-H]   Δ =   0.244 Å,  Initial =   0.859 Å

============================= Minor Angle Changes ==============================
Angle (5, 1, 63)   [N-Mn-H]  Δ =  16.471 °,  Initial =  96.799 °
Angle (61, 1, 63)  [C-Mn-H]  Δ =  15.528 °,  Initial =  81.202 °
Angle (2, 1, 63)   [P-Mn-H]  Δ =  13.032 °,  Initial = 171.266 °

Note: These angles depend on other changes and may not be significant alone.

============================ Minor Dihedral Changes ============================
Dihedral (63, 1, 2, 36)  [H-Mn-P-C]  Δ =  81.780 °,  Initial =  76.752 °

Note: These dihedrals depend on other changes and may not be significant alone.

================================================================================
```

**Interpretation:** Hydrogenation mechanism involving multiple N-H, H-O, and Mn-H bond changes. Note the handling of metal-ligand interactions and the lower magnitude H-H bond detection due to a secondary H threshold.

---

## Command Line Interface

### Basic Usage

```bash
vib_analysis <input_file> [options]
```

### Options
```text
> vib_analysis -h

usage: vib_analysis [-h] [--mode MODE] [--ts-frame TS_FRAME] [--bond-tolerance BOND_TOLERANCE] [--bond-threshold BOND_THRESHOLD] [--angle-threshold ANGLE_THRESHOLD]
                    [--dihedral-threshold DIHEDRAL_THRESHOLD] [--coupled-motion-filter COUPLED_MOTION_FILTER] [--all] [--graph] [--method {cheminf,xtb}] [--charge CHARGE]
                    [--multiplicity MULTIPLICITY] [--distance-tolerance DISTANCE_TOLERANCE] [--ascii-scale ASCII_SCALE] [--show-h] [--ascii-shells ASCII_SHELLS]
                    [--save-displacement] [--displacement-scale DISPLACEMENT_SCALE] [--no-save] [--orca-path ORCA_PATH] [--debug]
                    input

Analyze vibrational trajectories for structural changes

positional arguments:
  input                 Input file (XYZ trajectory or QM output)

options:
  -h, --help            show this help message and exit
  --mode, -m MODE       Vibrational mode to analyze (default: 0, ignored for XYZ)
  --ts-frame TS_FRAME   Frame index to use as TS reference (default: 0)
  --debug, -d           Enable debug output

vibrational analysis parameters:
  --bond-tolerance BOND_TOLERANCE
                        Bond detection tolerance factor (default: 1.4)
  --bond-threshold BOND_THRESHOLD
                        Threshold for significant bond changes in Å (default: 0.4)
  --angle-threshold ANGLE_THRESHOLD
                        Threshold for significant angle changes in degrees (default: 10.0)
  --dihedral-threshold DIHEDRAL_THRESHOLD
                        Threshold for significant dihedral changes in degrees (default: 20.0)
  --coupled-motion-filter COUPLED_MOTION_FILTER
                        Bond stability threshold for filtering coupled changes in Å (default: 0.2, advanced)
  --coupled-proton-threshold COUPLED_PROTON_THRESHOLD
                        Reduced threshold for coupled proton transfers in Å (default: 0.15, use "false" to disable)
  --all, -a             Report all changes including minor ones

graph analysis parameters:
  --graph, -g           Enable graph-based analysis and mode characterization (rotations, inversions, aromatic systems)
  --method {cheminf,xtb}
                        Graph building method (default: cheminf)
  --charge CHARGE       Molecular charge for graph building (default: 0)
  --multiplicity MULTIPLICITY
                        Spin multiplicity (auto-detected if not specified)
  --distance-tolerance DISTANCE_TOLERANCE
                        Tolerance for bond formation/breaking (default: 0.2 Å)
  --independent-graphs, -ig
                        Build molecular graphs from the displaced geometries rather than TS geometry with guided bonding (more rigorous for use with IRC or QRC displaced trajectories)

ASCII rendering options:
  --ascii-scale, -as ASCII_SCALE
                        Scale for ASCII molecular rendering (default: 2.5)
  --show-h              Show hydrogen atoms in ASCII rendering
  --ascii-shells, -ash ASCII_SHELLS
                        Neighbor shells around transformation core (default: 1)

output options:
  --save-displacement, -sd
                        Save displaced structure pair
  --displacement-scale, -ds DISPLACEMENT_SCALE
                        Displacement level (1-4, ~0.2-0.8 amplitude) (default: 1)
  --no-save             Do not save trajectory to disk (keep in memory only)
  --orca-path ORCA_PATH
                        Path to ORCA executable directory
```


### Threshold Tuning

```bash
# Adjust bond detection sensitivity
vib_analysis input.xyz --bond-threshold 0.3

# Adjust angle detection
vib_analysis input.xyz --angle-threshold 15.0
```

### Graph Analysis Options

```bash
# With ASCII visualization
vib_analysis input.xyz -g --ascii-scale 2.5 --show-h

# Adjust display around reactive center
vib_analysis input.xyz -g --ascii-shells 2

# Set molecular charge
vib_analysis input.xyz -g --charge -1

# Use independent graph building (more rigorous for IRC/QRC trajectories)
vib_analysis input.xyz -g --independent-graphs
```

### Independent Graph Building

By default, molecular graphs are built from TS geometry with bonding guided by the bond changes across the trajectory. The `--independent-graphs` flag enables a more rigorous approach where each displaced structure's graph is built independently using its actual geometry:

```bash
# Standard approach (TS-centric, default)
vib_analysis irc_trajectory.xyz -g

# Independent approach (builds from actual geometries)
vib_analysis irc_trajectory.xyz -g --independent-graphs
```

**When to use `--independent-graphs`:**
- Analyzing IRC or QRC trajectories with actual minima geometries
- Formal validation of connectivity changes
- Comparing results between methods
- Ensuring no bias from TS geometry

**Differences:**
- **Default (TS-centric)**: All graphs built from TS geometry with bond/unbond guidance based on distance changes
- **Independent**: Each graph built from its own actual geometry; `xyzgraph` determines connectivity from coordinates alone

### Output Control

```bash
# Save displaced structures
vib_analysis input.xyz --save-displacement --displacement-scale 2
# or
vib_analysis input.xyz -sd -ds 2


# Don't save trajectory to disk
vib_analysis input.xyz --no-save

# Specify ORCA path
vib_analysis input.out --orca-path /opt/orca
```

### Complete Example

```bash
vib_analysis bimp.out \
  --mode 0 \
  --graph \
  --debug \
  --save-displacement \
  --ascii-shells 1 \
  --ascii-scale 2.5
```

---

## Python API

See examples/examples.ipynb
This function will return a dictionary of the results, and printing can be turned on to produce the same as the CLI
For example:
```python
from vib_analysis import run_vib_analysis

orca_out = 'data/bimp.v000.xyz'

# Basic analysis
results = run_vib_analysis(
        input_file=orca_out,
    )

vib = results['vibrational']
print(vib)

theoretical_bond_changes = [(11,12), (10,14)]
if all(bond in vib['bond_changes'] for bond in theoretical_bond_changes):
    print(f'True: All theoretical bond changes {theoretical_bond_changes} found in results.')

# With graph analysis and independent graphs (for IRC/QRC trajectories)
results_ig = run_vib_analysis(
        input_file='irc_trajectory.xyz',
        enable_graph=True,
        independent_graphs=True,  # Build from actual geometries
        print_output=True
    )
```
Outputs:
```python
{'bond_changes': {(10, 14): (0.426, 2.656), (11, 12): (2.052, 2.064)}, 'angle_changes': {}, 'minor_angle_changes': {(13, 12, 29): (14.436, 122.116), (29, 12, 30): (12.54, 117.79), (12, 13, 14): (14.118, 123.702)}, 'dihedral_changes': {}, 'minor_dihedral_changes': {(0, 1, 10, 11): (36.48, 14.986), (4, 9, 10, 11): (50.966, 169.776), (29, 12, 13, 31): (67.358, 17.521), (12, 13, 31, 33): (62.151, 330.369)}, 'frame_indices': [5, 15], 'atom_index_map': {0: 'O', 1: 'C', ...}}
True: All theoretical bond changes [(11, 12), (10, 14)] found in results.
```
  - This can be used to check for a known vibrational mode (theoretical_bond_change) in `results['bond_changes']`
  - So in theory this could identify whether the correct TS mode has been identidied in a high throughput search if the atom indices are known (or available automatically)


### Results Structure

```python
{
    'trajectory': {
        'frames': List[Dict],        # List of frame dictionaries
        'frequencies': List[float],  # cm⁻¹ (None for XYZ)
        'trajectory_file': str       # Path to saved file
    },
    'vibrational': {
        'bond_changes': Dict[Tuple, Tuple[float, float]],
        'angle_changes': Dict[Tuple, Tuple[float, float]],
        'dihedral_changes': Dict[Tuple, Tuple[float, float]],
        'minor_angle_changes': Dict,
        'minor_dihedral_changes': Dict,
        'frame_indices': List[int],
        'atom_index_map': Dict[int, str]
    },
    'graph': {                       # Only if enable_graph=True
        'comparison': Dict,
        'ts_graph': nx.Graph,
        'frame1_graph': nx.Graph,
        'frame2_graph': nx.Graph,
        'ascii_ts': str,
        'ascii_ref': str,
        'ascii_disp': str
    },
    'displacement_files': Tuple[str, str]  # If save_displacement=True
}
```

---

## Advanced Options

### Configuration Parameters

All defaults are in `config.py` and can be overridden:

**Detection Tolerances:**
```python
BOND_TOLERANCE = 1.4        # vdW radii multiplier for TS
```

**Significance Thresholds:**
```python
BOND_THRESHOLD = 0.4        # Minimum Δ (Å) 
ANGLE_THRESHOLD = 10.0      # Minimum Δ (degrees)
DIHEDRAL_THRESHOLD = 20.0   # Minimum Δ (degrees)
COUPLED_MOTION_FILTER = 0.2 # For secondary filtering of coupled changes
COUPLED_PROTON_THRESHOLD = 0.15 # Low threshold for coupled H movements (Å) 
```

**Graph Analysis:**
```python
DISTANCE_TOLERANCE = 0.2    # Bond formation/breaking (Å)
ASCII_SCALE = 2.5           # Rendering scale
ASCII_NEIGHBOR_SHELLS = 1   # Expansion around reactive center
```

### Coupled Proton Transfer Detection

For systems involving proton transfers or H₂ coordination, a reduced threshold can detect coupled H movements that fall below the standard bond threshold allowing for chemically relevant but asynchronous bond changes:

**How it works:**
- When an H atom is involved in a detected bond change, all other bonds involving that H are checked with the reduced threshold (0.15 Å instead of 0.4 Å)

**CLI Usage:**
```bash
# Custom threshold (enabled by default)
vib_analysis input.xyz --coupled-proton-threshold 0.20

# Disable feature
vib_analysis input.xyz --coupled-proton-threshold false
```

**Python API:**
```python
# Custom threshold
results = run_vib_analysis('input.xyz', coupled_proton_threshold=0.20)

# Disabled
results = run_vib_analysis('input.xyz', coupled_proton_threshold=False)
```

### Displaced Structure Export

Generate structures for tight optimization to either side of the TS:

```bash
# Default: ±1 amplitude (~0.2)
vib_analysis input.xyz --save-displacement

# Higher amplitude: ±2 (~0.4)
vib_analysis input.xyz --save-displacement --displacement-scale 2

# Creates: input_F.xyz (forward), input_R.xyz (reverse)
```

Displacement scale 1-4 correspond to amplitudes of ~0.2, 0.4, 0.6, 0.8. The direction is arbitrary. 

### Custom Frame Selection

```bash
# Override TS frame
vib_analysis input.xyz --ts-frame 5
```

By default, frame 1 is the TS, and frames with maximum RMSD are selected automatically.

---

### Threshold Validation

The default bond displacement threshold (0.4 Å) has been validated against a diverse set of 15 transition state systems. "Ground truth" bond changes were determined from IRC calculations.

**Validation script:**
```bash
python examples/threshold_tuning.py
```

**Results:** The default threshold of 0.4 Å combined with coupled proton detection provides optimal performance:
- **100%** F1 score (balance of precision and recall)
- **100%** detection rate (all expected bonds found)
- **0%** false positive rate
- **full accuracy of *all* vibrational bonds across *all* 15 transition states**  

![threshold tuning](images/threshold_optimization.png)

Detailed validation results are written to `examples/threshold_optimization.txt` for full transparency.

---

## Important Notes

**Integrated threshold adjustments**
- Thresholds are reduced by 50% if there is no initial detection of interal coordinate changes
- This is flagged in the output and these may be less reliable
- Allows for the detection of changed coordinated in very low magnitude modes (*i.e.* hindered aryl rotation)

**Use as indicators only!** Always cross-validate with:
- IRC
- Optimisations of displaced structures
- Your own chemical insight

### File Formats

**Supported Inputs:**
- XYZ trajectory (`.xyz`) - direct read
- ORCA output (`.out`) - via cclib or orca_pltvib
- Gaussian (`.log`) - via cclib

**XYZ Format:**
```
<n_atoms>
Comment line
<symbol> <x> <y> <z>
...
```

Must contain ≥2 frames.

---

## Acknowledgments

- Uses [xyzgraph](https://github.com/aligfellow/xyzgraph) for graph construction and ascii printing
- QM output parsing via [cclib](https://github.com/cclib/cclib)
- Visualization examples with [v.2.0](https://github.com/briling/v) by Ksenia Briling

---
