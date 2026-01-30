"""
Utility functions for vibrational analysis.

This module provides file I/O utilities for trajectories and displaced structures,
geometry calculations, and centralized logging configuration.
"""

import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger("graphrc")

# Geometry calculation functions


def calculate_distance(positions: np.ndarray, i: int, j: int) -> float:
    """
    Calculate distance between atoms i and j.

    Args:
        positions: Nx3 array of atomic positions
        i, j: Atom indices

    Returns
    -------
        Distance in Angstroms
    """
    return round(float(np.linalg.norm(positions[j] - positions[i])), 3)


def calculate_angle(positions: np.ndarray, i: int, j: int, k: int) -> float:
    """
    Calculate angle between atoms i-j-k.

    Args:
        positions: Nx3 array of atomic positions
        i, j, k: Atom indices

    Returns
    -------
        Angle in degrees
    """
    v1 = positions[i] - positions[j]
    v2 = positions[k] - positions[j]

    cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_angle = np.clip(cos_angle, -1.0, 1.0)
    angle = np.arccos(cos_angle) * 180.0 / np.pi

    return round(float(angle), 3)


def calculate_dihedral(positions: np.ndarray, i: int, j: int, k: int, m: int) -> float:
    """
    Calculate dihedral angle for atoms i-j-k-m.

    Args:
        positions: Nx3 array of atomic positions
        i, j, k, m: Atom indices (dihedral convention i-j-k-l, using m for linting)

    Returns
    -------
        Dihedral angle in degrees
    """
    b1 = positions[j] - positions[i]
    b2 = positions[k] - positions[j]
    b3 = positions[m] - positions[k]

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    m1 = np.cross(n1, b2 / np.linalg.norm(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    dihedral = np.arctan2(y, x) * 180.0 / np.pi

    return round(float(dihedral), 3)


def write_xyz(filename: str, symbols: List[str], positions: np.ndarray, comment: str = "") -> None:
    """
    Write single structure to XYZ file.

    Args:
        filename: Output file path
        symbols: List of element symbols
        positions: Nx3 array of atomic positions
        comment: Comment line for XYZ file
    """
    with open(filename, "w") as f:
        f.write(f"{len(symbols)}\n")
        f.write(f"{comment}\n")
        for sym, pos in zip(symbols, positions):
            f.write(f"{sym} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")


def setup_logging(debug: bool = False) -> None:
    """
    Configure logging for graphRC package.

    Provides centralized logging setup for both CLI and API usage.
    In debug mode, prints a simple message before debug logging begins.

    Args:
        debug: Enable DEBUG level logging (also prints initialization message)
    """
    # Determine log level
    if debug:
        level = logging.DEBUG
        print("Initiating debugging:")
    else:
        level = logging.WARNING

    # Configure logging
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s", force=True)


def write_trajectory_file(trajectory_string: str, output_path: str) -> str:
    """
    Write trajectory string to XYZ file.

    Args:
        trajectory_string: Complete XYZ trajectory as string
        output_path: Path where file will be written

    Returns
    -------
        The output path (same as input)
    """
    with open(output_path, "w") as f:
        f.write(trajectory_string)
    logger.info(f"Wrote trajectory to {output_path}")
    return output_path


def write_displaced_structures(
    frames: List[Dict[str, Any]],
    prefix: str,
    indices: Optional[List[int]] = None,
    ts_frame: int = 0,
    overwrite: bool = True,
) -> List[str]:
    """
    Write displaced structures (forward/reverse) as XYZ files.

    Args:
        frames: List of frame dicts with 'symbols' and 'positions' keys
        prefix: Output filename prefix
        indices: Optional list of frame indices [forward, reverse]
        ts_frame: TS frame index (used if indices is None)
        overwrite: Whether to overwrite existing files

    Behavior:
      - If indices provided: use first → {prefix}_F.xyz, second → {prefix}_R.xyz
      - If indices is None: use (ts_frame-1, ts_frame+1)
      - Negative indices allowed (Python style)
      - Single index → only _F written

    Returns
    -------
        List of written file paths
    """
    written = []
    if not frames or len(frames) < 2:
        return written

    n = len(frames)

    def norm(i: int) -> Optional[int]:
        return i % n if -n <= i < n else None

    if indices is None:
        a = norm(ts_frame - 1)
        b = norm(ts_frame + 1)
        if a is None or b is None or a == b:
            a, b = 0, n - 1
        indices = [a, b]
    else:
        # Keep only first two valid normalized indices
        normed = []
        for raw in indices:
            ni = norm(raw)
            if ni is not None:
                normed.append(ni)
            if len(normed) == 2:
                break
        indices = normed

    if not indices:
        return written

    # First index (_F)
    f_idx = indices[0]
    f_path = f"{prefix}_F.xyz"
    if overwrite or not os.path.exists(f_path):
        frame = frames[f_idx]
        write_xyz(f_path, frame["symbols"], frame["positions"])
    written.append(f_path)

    # Optional second index (_R)
    if len(indices) > 1:
        r_idx = indices[1]
        if r_idx != f_idx:
            r_path = f"{prefix}_R.xyz"
            if overwrite or not os.path.exists(r_path):
                frame = frames[r_idx]
                write_xyz(r_path, frame["symbols"], frame["positions"])
            written.append(r_path)

    return written


def save_displacement_pair(
    frames: List[Dict[str, Any]],
    ts_frame: int,
    output_prefix: str,
    scale: int = 1,
    max_level: int = 4,
    print_output: bool = False,
) -> Optional[Tuple[str, str]]:
    """
    Save symmetric displaced structures (TS ± scale) as XYZ files.

    Uses Python-style negative index wrapping (e.g., frame -1 is the last frame).

    Args:
        frames: List of frame dicts with 'symbols' and 'positions' keys
        ts_frame: Index of TS frame
        output_prefix: Prefix for output files
        scale: Displacement scale (1-4, corresponding to amplitudes ~0.2-0.8)
        max_level: Maximum allowed scale
        print_output: Print status messages

    Returns
    -------
        (forward_path, reverse_path) if successful, else None
    """
    n_frames = len(frames)

    if not (1 <= scale <= max_level):
        logger.warning(f"Invalid scale {scale} (must be 1-{max_level})")
        if print_output:
            print(f"Invalid scale {scale} (must be 1-{max_level}).")
        return None

    # Calculate indices with wrapping support
    minus_idx = ts_frame - scale
    plus_idx = ts_frame + scale

    # Normalize indices (allow Python-style negative indexing)
    def normalize_index(idx: int) -> Optional[int]:
        """Normalize index to valid range, supporting negative indices."""
        if -n_frames <= idx < n_frames:
            return idx % n_frames
        return None

    norm_minus = normalize_index(minus_idx)
    norm_plus = normalize_index(plus_idx)

    if norm_minus is None or norm_plus is None:
        logger.warning(f"Scale {scale} out of range for TS {ts_frame} (total {n_frames} frames)")
        if print_output:
            print(f"Scale {scale} out of range for TS {ts_frame} (total {n_frames} frames).")
        return None

    paths = write_displaced_structures(frames, prefix=output_prefix, indices=[norm_minus, norm_plus])

    if len(paths) == 2:
        logger.info(f"Saved displaced pair (±{scale}): {os.path.basename(paths[0])}, {os.path.basename(paths[1])}")
        if print_output:
            print(f"Saved displaced pair (±{scale}): {os.path.basename(paths[0])}, {os.path.basename(paths[1])}")
        return (paths[0], paths[1])

    return None
