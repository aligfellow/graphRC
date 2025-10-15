"""
Utility functions for vibrational analysis.

This module provides file I/O utilities for trajectories and displaced structures,
as well as centralized logging configuration.
"""

import os
import logging
from typing import List, Optional, Tuple
from ase import Atoms
from ase.io import write

logger = logging.getLogger("vib_analysis")


def setup_logging(verbose: bool = False, debug: bool = False) -> None:
    """
    Configure logging for vib_analysis package.
    
    Provides centralized logging setup for both CLI and API usage.
    In debug mode, prints a simple message before debug logging begins.
    
    Args:
        verbose: Enable INFO level logging
        debug: Enable DEBUG level logging (also prints initialization message)
    """
    # Determine log level
    if debug:
        level = logging.DEBUG
        print("Initiating debugging:")
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    
    # Configure logging
    logging.basicConfig(
        level=level,
        format='[%(levelname)s] %(message)s',
        force=True 
    )


def write_trajectory_file(trajectory_string: str, output_path: str) -> str:
    """
    Write trajectory string to XYZ file.
    
    Args:
        trajectory_string: Complete XYZ trajectory as string
        output_path: Path where file will be written
        
    Returns:
        The output path (same as input)
    """
    with open(output_path, 'w') as f:
        f.write(trajectory_string)
    logger.info(f"Wrote trajectory to {output_path}")
    return output_path


def write_displaced_structures(
    frames: List[Atoms],
    prefix: str,
    indices: Optional[List[int]] = None,
    ts_frame: int = 0,
    overwrite: bool = True
) -> List[str]:
    """
    Write displaced structures (forward/reverse) as XYZ files.
    
    Args:
        frames: List of ASE Atoms objects
        prefix: Output filename prefix
        indices: Optional list of frame indices [forward, reverse]
        ts_frame: TS frame index (used if indices is None)
        overwrite: Whether to overwrite existing files
    
    Behavior:
      - If indices provided: use first → {prefix}_F.xyz, second → {prefix}_R.xyz
      - If indices is None: use (ts_frame-1, ts_frame+1)
      - Negative indices allowed (Python style)
      - Single index → only _F written
    
    Returns:
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
        write(f_path, frames[f_idx], format='xyz')
    written.append(f_path)
    
    # Optional second index (_R)
    if len(indices) > 1:
        r_idx = indices[1]
        if r_idx != f_idx:
            r_path = f"{prefix}_R.xyz"
            if overwrite or not os.path.exists(r_path):
                write(r_path, frames[r_idx], format='xyz')
            written.append(r_path)
    
    return written


def save_displacement_pair(
    frames: List[Atoms],
    ts_frame: int,
    output_prefix: str,
    scale: int = 1,
    max_level: int = 4,
    verbose: bool = False,
) -> Optional[Tuple[str, str]]:
    """
    Save symmetric displaced structures (TS ± level) as XYZ files.
    
    Args:
        frames: List of trajectory frames
        ts_frame: Index of TS frame
        output_prefix: Prefix for output files
        level: Displacement level (1-4, corresponding to amplitudes ~0.2-0.8)
        max_level: Maximum allowed level
        verbose: Print status messages
    
    Returns:
        (forward_path, reverse_path) if successful, else None
    """
    n_frames = len(frames)
    
    if not (1 <= level <= max_level):
        logger.warning(f"Invalid level {level} (must be 1–{max_level})")
        if verbose:
            print(f"Invalid level {level} (must be 1–{max_level}).")
        return None
    
    minus_idx = ts_frame - level
    plus_idx = ts_frame + level
    
    if not (0 <= minus_idx < n_frames and 0 <= plus_idx < n_frames):
        logger.warning(
            f"Level {level} out of range for TS {ts_frame} "
            f"(total {n_frames} frames)"
        )
        if verbose:
            print(
                f"Level {level} out of range for TS {ts_frame} "
                f"(total {n_frames} frames)."
            )
        return None
    
    paths = write_displaced_structures(
        frames, 
        prefix=output_prefix, 
        indices=[minus_idx, plus_idx]
    )
    
    if len(paths) == 2:
        logger.info(
            f"Saved displaced pair (±{level}): "
            f"{os.path.basename(paths[0])}, {os.path.basename(paths[1])}"
        )
        if verbose:
            print(
                f"Saved displaced pair (±{level}): "
                f"{os.path.basename(paths[0])}, {os.path.basename(paths[1])}"
            )
    
    return tuple(paths) if len(paths) == 2 else None
