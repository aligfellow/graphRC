from importlib.metadata import version

__version__ = version("graphRC")
__citation__ = f"A. S. Goodfellow, graphRC: Internal Coordinate Analysis of Vibrational Modes, v{__version__}, 2025, https://github.com/aligfellow/graphRC.git."

from .api import load_trajectory, run_vib_analysis
