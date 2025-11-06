from importlib.metadata import version

__version__ = version("vib_analysis")
__citation__ = f"A. S. Goodfellow, vib_analysis: Internal Coordinate Analysis of Vibrational Modes, v{__version__}, 2025, https://github.com/aligfellow/vib_analysis.git."

from .api import run_vib_analysis, load_trajectory
