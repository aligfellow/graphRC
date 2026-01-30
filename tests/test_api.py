"""Smoke tests for graphrc API."""

from pathlib import Path

from graphrc import __citation__, __version__

DATA_DIR = Path(__file__).parent.parent / "examples" / "data"
TEST_FILE = str(DATA_DIR / "sn2.v000.xyz")


def test_version():
    """Version is defined."""
    assert __version__


def test_citation():
    """Citation contains package name."""
    assert "graphRC" in __citation__

