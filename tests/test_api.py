"""Smoke tests for graphrc API."""

import shutil
from pathlib import Path
from unittest.mock import Mock

import numpy as np
import pytest

from graphrc import __citation__, __version__, analyse, displace, run_vib_analysis, trajectory
from graphrc.convert import (
    _make_amplitudes,
    parse_cclib_output,
    parse_orca_output,
    parse_xyz_string_to_frames,
    validate_n_frames,
)
from graphrc.output import print_analysis_results, print_metadata_header
from graphrc.utils import save_displacement_pair

DATA_DIR = Path(__file__).resolve().parents[1] / "examples" / "data"


def test_version():
    """Version is defined."""
    assert __version__


def test_citation():
    """Citation contains package name."""
    assert "graphRC" in __citation__


def test_make_amplitudes_default():
    """The 20-frame default reproduces the original hardcoded amplitude list."""
    expected = [
        0.0, -0.2, -0.4, -0.6, -0.8, -1.0, -0.8, -0.6, -0.4, -0.2,
        0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.8, 0.6, 0.4, 0.2,
    ]  # fmt: skip
    assert np.allclose(_make_amplitudes(20), expected)
    assert len(_make_amplitudes(40)) == 40


def test_validate_n_frames():
    """n_frames must be a positive multiple of 4."""
    assert validate_n_frames(20) == 20
    for bad in (0, 7, -4):
        with pytest.raises(ValueError, match="multiple of 4"):
            validate_n_frames(bad)


@pytest.fixture(
    scope="module",
    params=[("sn2.out", parse_orca_output), ("spiro-ts1.log", parse_cclib_output), ("dihedral.out", parse_orca_output)],
    ids=["orca-reaction", "gaussian-reaction", "orca-rotation"],
)
def parsed_mode(request):
    name, parser = request.param
    path = DATA_DIR / name
    frequencies, xyz, mode = parser(str(path), 0)
    return path, frequencies, parse_xyz_string_to_frames(xyz), mode


def test_trajectory_matches_parsed_mode(parsed_mode):
    _, _, parsed_frames, mode = parsed_mode
    equilibrium = parsed_frames[0]
    frames = trajectory(equilibrium["symbols"], equilibrium["positions"], mode)
    assert len(frames) == len(parsed_frames)
    for actual, expected in zip(frames, parsed_frames):
        assert actual["symbols"] == expected["symbols"]
        np.testing.assert_allclose(actual["positions"], expected["positions"], rtol=0, atol=5.1e-7)


def test_analyse_matches_file_pipeline(parsed_mode, tmp_path, monkeypatch):
    path, frequencies, parsed_frames, mode = parsed_mode
    monkeypatch.chdir(tmp_path)
    expected = run_vib_analysis(str(path), enable_graph=True, save_trajectory=False)
    equilibrium = parsed_frames[0]
    actual = analyse(
        equilibrium["symbols"], equilibrium["positions"], mode, frequency=frequencies[0], enable_graph=True
    )
    assert actual.keys() == expected.keys()
    assert actual["vibrational"] == expected["vibrational"]
    assert actual["characterization"] == expected["characterization"]
    assert actual["graph"]["comparison"] == expected["graph"]["comparison"]
    assert actual["trajectory"]["frequencies"] == [frequencies[0]]
    assert actual["trajectory"]["trajectory_file"] is None
    assert actual["metadata"]["input_file"] is None
    assert list(tmp_path.iterdir()) == []


def test_analyse_forwards_distance_tolerance(monkeypatch):
    graph_analysis = Mock(return_value={})
    monkeypatch.setattr("graphrc.api.analyze_displacement_graphs", graph_analysis)
    analyse(
        ["H", "H"],
        [[0, 0, 0], [0, 0, 0.74]],
        [[0, 0, 0], [0, 0, 0.5]],
        enable_graph=True,
        distance_tolerance=0.15,
    )
    assert graph_analysis.call_args.kwargs["distance_tolerance"] == 0.15


@pytest.mark.parametrize("frequency", [None, 0.0, -100.0])
def test_array_api_preserves_inputs_and_writes_no_files(frequency, tmp_path, monkeypatch, capsys):
    monkeypatch.chdir(tmp_path)
    positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.740000123]])
    mode = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.5]])
    original_positions, original_mode = positions.copy(), mode.copy()
    frames = trajectory([1, 1], positions, mode, n_frames=4)
    np.testing.assert_array_equal(frames[0]["positions"], positions)
    np.testing.assert_array_equal(frames[1]["positions"], positions - mode)
    np.testing.assert_array_equal(frames[3]["positions"], positions + mode)
    assert frames[0]["symbols"] == ["H", "H"]
    displace(positions, mode)
    results = analyse([1, 1], positions, mode, frequency=frequency)
    assert results["vibrational"]["bond_changes"] == {(0, 1): (1.0, 0.74)}
    assert results["trajectory"]["frequencies"] == (None if frequency is None else [frequency])
    assert results["characterization"] is None
    assert results["graph"] is None
    assert results["displacement_files"] is None
    assert capsys.readouterr().out == ""
    print_metadata_header(results["metadata"], results["trajectory"])
    print_analysis_results(results)
    assert "Input:          in memory" in capsys.readouterr().out
    assert list(tmp_path.iterdir()) == []
    frames[0]["positions"][:] = 99
    results["trajectory"]["mode_displacement"][:] = 99
    np.testing.assert_array_equal(positions, original_positions)
    np.testing.assert_array_equal(mode, original_mode)


@pytest.mark.parametrize("scale", [1, 2, 6])
def test_displace_matches_saved_pair(scale, tmp_path):
    positions = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.74]])
    mode = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]])
    forward, reverse = displace(positions, mode, scale=scale)
    paths = save_displacement_pair(
        trajectory(["H", "H"], positions, mode),
        ts_frame=0,
        output_prefix=str(tmp_path / "pair"),
        scale=scale,
        mode_displacement=mode,
    )
    assert paths is not None
    for actual, path in zip((forward, reverse), paths):
        saved = parse_xyz_string_to_frames(Path(path).read_text())[0]
        np.testing.assert_allclose(actual, saved["positions"], rtol=0, atol=5.1e-7)


@pytest.mark.parametrize(
    ("positions", "mode", "message"),
    [
        ([0, 0, 0], [0, 0, 0], "positions must have shape"),
        (np.zeros((2, 2)), np.zeros((2, 2)), "positions must have shape"),
        (np.empty((0, 3)), np.empty((0, 3)), "at least one atom"),
        (np.zeros((2, 3)), np.zeros((1, 3)), "does not match"),
        (np.zeros((2, 3)), np.zeros(6), "does not match"),
        ([[np.nan, 0, 0]], [[0, 0, 0]], "finite values"),
        ([[0, 0, 0]], [[np.inf, 0, 0]], "finite values"),
    ],
)
def test_invalid_mode_arrays(positions, mode, message):
    with pytest.raises(ValueError, match=message):
        displace(positions, mode)


def test_trajectory_validates_symbols_and_frame_count():
    positions, mode = np.zeros((2, 3)), np.zeros((2, 3))
    with pytest.raises(ValueError, match="does not match"):
        analyse(["H", "H"], positions, np.zeros((1, 3)))
    for symbols, message in [
        (["H"], "one element per atom"),
        (["H", "bad"], "valid element"),
        ([1, 999], "valid element"),
    ]:
        with pytest.raises(ValueError, match=message):
            trajectory(symbols, positions, mode)
    with pytest.raises(ValueError, match="multiple of 4"):
        trajectory(["H", "H"], positions, mode, n_frames=7)
    with pytest.raises(ValueError, match="finite scalar"):
        displace(positions, mode, scale=np.inf)


def test_file_api_keeps_default_trajectory_output(tmp_path):
    path = tmp_path / "sn2.out"
    shutil.copyfile(DATA_DIR / "sn2.out", path)
    results = run_vib_analysis(str(path))
    output = path.with_suffix(".v000.xyz")
    assert results["trajectory"]["trajectory_file"] == str(output)
    assert output.is_file()
    assert results["displacement_files"] is None
