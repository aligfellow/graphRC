"""
Ground truth definitions for transition state validation.

This file contains the expected bond, angle, and dihedral changes
for each test system in examples/data/

Atom indices are zero-indexed.
"""

EXPECTED_RESULTS = {
    'ambh': {
        'bonds': [(2, 11), (11, 18)] # C-H and H-N
    },
    'atrop': {
        'dihedrals': [(1, 5, 6, 7), (6, 7, 25, 27)]
    },
    'bimp': {
        'bonds': [(11, 12), (10, 14)],  # O-C and C-C
    },
    'cpa': {
        'bonds': [(23, 35), (35, 75), (70, 74)] # O-H, H-O, C-C
    },
    'cycl': {
        'bonds': [(2, 6)],  # C-N cyclisation
    },
    'dihedral': {
        'dihedrals': [(6, 0, 3, 7)],  # F-C-C-F rotation
    },
    'mn-h2': {
        'bonds': [(1, 64), (5, 65), (63, 64), (64, 66), (65, 66)],  # Mn-H, N-H, H-H, H-O, H-O
    },
    'mn-hy': {
        'bonds': [(1, 67), (67, 77)] # Mn-H, H-C
    },
    'ru-co-nuc': {
        'bonds': [(3, 66), (64, 65), (65, 66)],  # N-H, C-N, N-H
    },
    'ru-co': {
        'bonds': [(0, 63), (0, 64), (3, 66), (64, 66)],  # Ru-O, Ru-C, N-H, C-H
    },
    'sn2': {
        'bonds': [(0, 4), (0, 5)],  # C-F and C-Cl
    },
    'sn2_large': {
        'bonds': [(0, 1), (0, 21)],  # C-I and C-N
    },
    'spiro-ts1': {
        'bonds': [(14, 15)],  # C-N
    },
    'spiro-ts2': {
        'bonds': [(13, 14)],  # O-C
    },
    'thia-ma': {
        'bonds': [(35, 47)],  # C-S
    },
}
