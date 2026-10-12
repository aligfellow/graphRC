"""
Trajectory conversion utilities for QM output files.

Handles ORCA and cclib-compatible formats.
"""

import logging
from typing import Any, Dict, List

import numpy as np
from xyzgraph import DATA

from . import config

logger = logging.getLogger("graphrc")


def validate_n_frames(n_frames: int) -> int:
    """Return *n_frames* unchanged if it is a positive multiple of 4, else raise ValueError.

    The trajectory amplitudes form a triangle wave; a multiple of 4 ensures the
    extremes land exactly on ±1.0 and the midpoint lands exactly on 0.0.
    """
    if n_frames < 4 or n_frames % 4 != 0:
        raise ValueError(f"n_frames must be a positive multiple of 4, got {n_frames}.")
    return n_frames


def _make_amplitudes(n_frames: int) -> np.ndarray:
    """Generate *n_frames* amplitude samples covering one full oscillation cycle.

    The cycle runs 0 → -1 → 0 → +1 → (wraps back to 0), sampled as a
    triangle wave.  n_frames must be a multiple of 4 so that the extremes
    land exactly on ±1.0 and the midpoint lands exactly on 0.0.
    """
    validate_n_frames(n_frames)
    t = np.linspace(0, 1, n_frames, endpoint=False)
    return np.where(t < 0.25, -4 * t, np.where(t < 0.75, 4 * t - 2, 4 - 4 * t))


def _validate_mode(positions, mode):
    """Validate Cartesian coordinates and mode without allowing broadcasting."""
    positions = np.asarray(positions, dtype=float)
    mode = np.asarray(mode, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 3 or len(positions) == 0:
        raise ValueError("positions must have shape (n_atoms, 3), with at least one atom")
    if mode.shape != positions.shape:
        raise ValueError(f"mode shape {mode.shape} does not match positions shape {positions.shape}")
    if not np.isfinite(positions).all() or not np.isfinite(mode).all():
        raise ValueError("positions and mode must contain only finite values")
    return positions, mode


def trajectory(symbols, positions, mode, n_frames: int = config.DEFAULT_N_FRAMES) -> List[Dict[str, Any]]:
    """Return frame dictionaries from element symbols or atomic numbers, without file I/O.

    Positions and mode are (n_atoms, 3) arrays in Å, used without normalization.
    Frames span ``positions ± mode``; n_frames must be a positive multiple of 4.
    """
    amplitudes = _make_amplitudes(n_frames)
    positions, mode = _validate_mode(positions, mode)
    symbols = [DATA.n2s.get(symbol, symbol) for symbol in symbols]
    if len(symbols) != len(positions):
        raise ValueError("symbols must contain one element per atom")
    if any(symbol not in DATA.s2n for symbol in symbols):
        raise ValueError("symbols must contain valid element symbols or atomic numbers")
    return [{"symbols": symbols.copy(), "positions": positions + amp * mode} for amp in amplitudes]


def displace(positions, mode, scale: float = 1.0):
    """Return forward and reverse arrays using ``positions ± 0.2 * scale * mode``.

    Positions and mode must both have shape (n_atoms, 3), in Å. No mode
    normalization is applied: a mode whose largest atomic displacement is 1 Å
    gives a largest step of 0.2 Å at scale=1. Inputs are never modified.
    """
    positions, mode = _validate_mode(positions, mode)
    if not np.isscalar(scale) or not np.isfinite(scale):
        raise ValueError("scale must be a finite scalar")
    step = 0.2 * scale * mode
    return positions + step, positions - step


def _trajectory_to_xyz(frames, mode, frequency):
    """Format generated frames with the existing QM trajectory comments and precision."""
    trj_data = ""
    for amp, frame in zip(_make_amplitudes(len(frames)), frames):
        trj_data += f"{len(frame['symbols'])}\n"
        trj_data += f"Mode: {mode}, Frequency: {frequency:.2f} cm**-1, Amplitude: {amp:.2f}\n"
        for sym, coord in zip(frame["symbols"], frame["positions"]):
            trj_data += f"{sym} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n"
    return trj_data


def is_orca_output(filepath: str) -> bool:
    """Return True if the file contains the ORCA banner in its header."""
    try:
        with open(filepath) as f:
            for _ in range(10):
                if "O   R   C   A" in f.readline():
                    return True
    except OSError:
        pass
    return False


def parse_cclib_output(output_file, mode, n_frames: int = config.DEFAULT_N_FRAMES):
    """
    Parse QM output with cclib and generate trajectory.

    Returns: (frequencies, trajectory_xyz_string).
    """
    mode = int(mode)
    validate_n_frames(n_frames)

    import cclib

    try:
        parser = cclib.io.ccopen(output_file)
        if parser is None:
            raise ValueError(f"cclib could not parse {output_file}.")
    except Exception as e:
        raise ValueError(f"Error parsing {output_file} with cclib: {e}") from e

    data = parser.parse()

    if not hasattr(data, "vibfreqs") or len(data.vibfreqs) == 0:
        raise ValueError("No vibrational frequencies found in file.")

    freqs = data.vibfreqs
    num_modes = len(freqs)
    if mode < 0 or mode >= num_modes:
        raise ValueError(f"Mode index {mode} out of range. File has {num_modes} modes.")

    atom_numbers = data.atomnos
    atom_symbols = [DATA.n2s[z] for z in atom_numbers]
    eq_coords = data.atomcoords[-1]
    displacement = np.array(data.vibdisps[mode])
    frames = trajectory(atom_symbols, eq_coords, displacement, n_frames)
    trj_data = _trajectory_to_xyz(frames, mode, freqs[mode])

    return freqs, trj_data, displacement


def parse_orca_output(orca_file: str, mode: int, n_frames: int = config.DEFAULT_N_FRAMES):
    """
    Parse ORCA output file directly for normal modes.

    Reads VIBRATIONAL FREQUENCIES, CARTESIAN COORDINATES (ANGSTROEM), and
    NORMAL MODES sections from the last frequency block in the file.

    The NORMAL MODES table stores mass-weighted Cartesian displacements
    (L_ik = disp_ik / sqrt(m_i)), normalised to unit length. We apply them
    with the same amplitude sequence as parse_cclib_output so that trajectories
    are directly comparable.

    Returns: (frequencies, trajectory_xyz_string).
    """
    mode = int(mode)
    validate_n_frames(n_frames)

    with open(orca_file) as f:
        lines = f.readlines()

    # ------------------------------------------------------------------ freqs
    freq_block_starts = [i for i, line in enumerate(lines) if "VIBRATIONAL FREQUENCIES" in line]
    if not freq_block_starts:
        raise ValueError("No VIBRATIONAL FREQUENCIES section found.")
    freq_start = freq_block_starts[-1]

    all_freqs: List[float] = []
    n_zero = 0
    for line in lines[freq_start:]:
        if "NORMAL MODES" in line:
            break
        parts = line.split(":")
        if len(parts) > 1:
            try:
                f = float(parts[1].split()[0])
                all_freqs.append(f)
                if abs(f) <= 1e-5:
                    n_zero += 1
            except (ValueError, IndexError):
                continue

    freqs = [f for f in all_freqs if abs(f) > 1e-5]
    if not freqs:
        raise ValueError("No non-zero vibrational frequencies found.")
    if mode < 0 or mode >= len(freqs):
        raise ValueError(f"Mode {mode} out of range — file has {len(freqs)} modes.")

    # column index in the NORMAL MODES table
    col = mode + n_zero

    # --------------------------------------------------------- equilibrium geometry
    # Use last CARTESIAN COORDINATES (ANGSTROEM) block before the last freq section
    coord_starts = [i for i, line in enumerate(lines) if "CARTESIAN COORDINATES (ANGSTROEM)" in line]
    coord_starts_before = [i for i in coord_starts if i < freq_start]
    if not coord_starts_before:
        coord_starts_before = coord_starts
    coord_start = coord_starts_before[-1] + 2  # skip header + dashes

    atom_symbols: List[str] = []
    eq_coords: List[List[float]] = []
    for line in lines[coord_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("-"):
            break
        parts = stripped.split()
        if len(parts) < 4:
            break
        try:
            atom_symbols.append(parts[0])
            eq_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        except ValueError:
            break

    if not atom_symbols:
        raise ValueError("Could not parse equilibrium geometry.")
    n_atoms = len(atom_symbols)
    n_dof = 3 * n_atoms
    eq = np.array(eq_coords)

    # -------------------------------------------------------------- NORMAL MODES
    nm_starts = [i for i, line in enumerate(lines) if line.strip() == "NORMAL MODES"]
    if not nm_starts:
        raise ValueError("No NORMAL MODES section found.")
    nm_start = nm_starts[-1]

    # skip header lines (dashes + description + blank)
    data_start = nm_start + 1
    while data_start < len(lines) and not lines[data_start].strip().startswith("0"):
        data_start += 1

    # parse table: blocks of 6 columns, n_dof rows each
    mode_col: List[float] = [0.0] * n_dof
    i = data_start
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        # header line: space-separated integers
        parts = line.split()
        try:
            col_indices = [int(x) for x in parts]
        except ValueError:
            break  # reached next section
        i += 1
        if col not in col_indices:
            # skip this block's rows
            i += n_dof
            continue
        local_idx = col_indices.index(col)
        for _ in range(n_dof):
            if i >= len(lines):
                break
            row_parts = lines[i].split()
            row_idx = int(row_parts[0])
            mode_col[row_idx] = float(row_parts[1 + local_idx])
            i += 1
        break  # found and parsed our column

    displacement = np.array(mode_col).reshape(n_atoms, 3)

    # --------------------------------------------------------- build trajectory
    frames = trajectory(atom_symbols, eq, displacement, n_frames)
    trj_data = _trajectory_to_xyz(frames, mode, freqs[mode])

    return freqs, trj_data, displacement


def parse_xyz_string_to_frames(trj_data_str: str) -> List[Dict[str, Any]]:
    """
    Parse XYZ trajectory string into list of frame dicts.

    Returns list of frames.
    """
    lines = trj_data_str.strip().split("\n")
    frames = []
    i = 0
    while i < len(lines):
        try:
            num_atoms = int(lines[i].strip())
        except (ValueError, IndexError):
            break
        if i + 1 >= len(lines):
            break
        lines[i + 1]
        coords = []
        symbols = []
        for j in range(num_atoms):
            if i + 2 + j >= len(lines):
                break
            parts = lines[i + 2 + j].split()
            symbols.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        frame = {"symbols": symbols, "positions": np.array(coords)}
        frames.append(frame)
        i += 2 + num_atoms
    return frames
