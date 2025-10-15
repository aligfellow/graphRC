"""
Trajectory conversion utilities for QM output files.
Handles ORCA and cclib-compatible formats.
"""
import os
import subprocess
import numpy as np
import logging
from ase import Atoms
from ase.io import write
import cclib
from typing import List, Tuple, Optional

logger = logging.getLogger("vib_analysis")


def get_orca_pltvib_path():
    """Find orca_pltvib executable in the same directory as orca"""
    orca_path = os.popen('which orca').read().strip()
    if not orca_path:
        raise RuntimeError("ORCA not found in PATH. Please ensure ORCA is installed.")
    
    orca_dir = os.path.dirname(orca_path)
    pltvib_path = os.path.join(orca_dir, 'orca_pltvib')
    
    if not os.path.exists(pltvib_path):
        raise RuntimeError(f"orca_pltvib not found at {pltvib_path}")
    
    return pltvib_path


def get_orca_frequencies(orca_file):
    """Extract vibrational frequencies from ORCA output"""
    with open(orca_file, 'r') as f:
        lines = f.readlines()
    
    section_indices = [i for i, line in enumerate(lines) if "VIBRATIONAL FREQUENCIES" in line]
    if not section_indices:
        raise ValueError("No vibrational frequencies section found in ORCA output.")
    
    idx = section_indices[-1]  # last occurrence if multiple are present
    freqs = []
    
    for line in lines[idx:]:
        if "NORMAL MODES" in line:
            break
        parts = line.split(':')
        if len(parts) > 1:
            try:
                freq = float(parts[1].split()[0])
                freqs.append(freq)
            except (ValueError, IndexError):
                continue
    freqs = [f for f in freqs if abs(f) > 1e-5]
    return freqs


def convert_orca(orca_file, mode, pltvib_path=None):
    """
    Convert ORCA output to vibration trajectory string.
    Returns: (frequencies, trajectory_xyz_string)
    """
    if pltvib_path is None:
        pltvib_path = get_orca_pltvib_path()
    if not os.path.exists(orca_file):
        raise FileNotFoundError(f"ORCA output file {orca_file} does not exist.")
    
    basename = os.path.splitext(orca_file)[0]
    
    # Determine orca_mode offset (5 or 6)
    with open(orca_file, 'r') as f:
        lines = f.readlines()
    
    n_atoms = None
    for line in lines:
        if "Number of atoms" in line:
            n_atoms = int(line.split()[-1])
            break
    if n_atoms is None:
        raise ValueError("Could not determine number of atoms from ORCA output.")
    
    orca_mode = int(mode) + (5 if n_atoms < 3 else 6)
    
    # Handle multiple frequency blocks
    freq_indices = [i for i, line in enumerate(lines) if "VIBRATIONAL FREQUENCIES" in line]
    if not freq_indices:
        raise ValueError("No vibrational frequencies section found in ORCA output.")
    
    coord_indices = [i for i, line in enumerate(lines) if "CARTESIAN COORDINATES (ANGSTROEM)" in line]
    
    if len(freq_indices) > 1:
        print(f"INFO: Multiple 'VIBRATIONAL FREQUENCIES' sections found. Using the last one.")
        idx = max(i for i in coord_indices if i < freq_indices[-1])
        tmp_file = f'{basename}.tmp'
        with open(tmp_file, 'w') as f:
            f.writelines(lines[idx:])
        subprocess.run([pltvib_path, tmp_file, str(orca_mode)], 
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        os.remove(tmp_file)
        os.system(f'mv {basename}.tmp.v{orca_mode:03d}.xyz {basename}.out.v{orca_mode:03d}.xyz')
    else:
        subprocess.run([pltvib_path, orca_file, str(orca_mode)], 
                      stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    
    os.system(f'mv {basename}.out.v{orca_mode:03d}.xyz {basename}.out.v{mode:03d}.xyz')
    orca_vib = f'{basename}.out.v{mode:03d}.xyz'
    
    if not os.path.exists(orca_vib):
        raise FileNotFoundError(f"File {orca_vib} not found. Ensure ORCA output is correct.")
    
    with open(orca_vib, 'r') as f:
        lines = f.readlines()
    
    xyz_len = int(lines[0].split()[0]) + 2
    xyzs = [lines[i:i+xyz_len] for i in range(0, len(lines), xyz_len)]
    
    trj_data = ""
    for idx_block in xyzs:
        trj_data += idx_block[0]
        trj_data += f"Mode {mode} Frame: {idx_block[1]}"
        for line in idx_block[2:]:
            parts = line.split()
            trj_data += f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}\n"
    
    os.remove(orca_vib)
    freqs = get_orca_frequencies(orca_file)
    return freqs, trj_data


def parse_cclib_output(output_file, mode):
    """
    Parse QM output with cclib and generate trajectory.
    Returns: (frequencies, trajectory_xyz_string)
    """
    mode = int(mode)
    amplitudes = [0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -0.8, -0.6, -0.4, -0.2, 
                  0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6, 0.4, 0.2]
    
    try:
        parser = cclib.io.ccopen(output_file)
        if parser is None:
            raise ValueError(f"cclib could not parse {output_file}.")
    except Exception as e:
        raise ValueError(f"Error parsing {output_file} with cclib: {e}")
    
    data = parser.parse()
    
    if not hasattr(data, 'vibfreqs') or len(data.vibfreqs) == 0:
        raise ValueError("No vibrational frequencies found in file.")
    
    freqs = data.vibfreqs
    num_modes = len(freqs)
    if mode < 0 or mode >= num_modes:
        raise ValueError(f"Mode index {mode} out of range. File has {num_modes} modes.")
    
    atom_numbers = data.atomnos
    atom_symbols = [Atoms(numbers=[z]).get_chemical_symbols()[0] for z in atom_numbers]
    eq_coords = data.atomcoords[-1]
    displacement = np.array(data.vibdisps[mode])
    freq = freqs[mode]
    
    trj_data = ""
    for amp in amplitudes:
        displaced = eq_coords + amp * displacement
        trj_data += f"{len(atom_numbers)}\n"
        trj_data += f"Mode: {mode}, Frequency: {freq:.2f} cm**-1, Amplitude: {amp:.2f}\n"
        for sym, coord in zip(atom_symbols, displaced):
            trj_data += f"{sym} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n"
    
    return freqs, trj_data


def parse_xyz_string_to_frames(trj_data_str: str) -> List[Atoms]:
    """
    Parse XYZ trajectory string into list of ASE Atoms objects.
    Returns list of frames.
    """
    lines = trj_data_str.strip().split('\n')
    frames = []
    i = 0
    while i < len(lines):
        try:
            num_atoms = int(lines[i].strip())
        except (ValueError, IndexError):
            break
        if i + 1 >= len(lines):
            break
        comment = lines[i + 1]
        coords = []
        symbols = []
        for j in range(num_atoms):
            if i + 2 + j >= len(lines):
                break
            parts = lines[i + 2 + j].split()
            symbols.append(parts[0])
            coords.append([float(x) for x in parts[1:4]])
        frame = Atoms(symbols=symbols, positions=coords)
        frames.append(frame)
        i += 2 + num_atoms
    return frames
