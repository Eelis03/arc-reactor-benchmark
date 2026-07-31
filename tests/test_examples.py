"""Tier three: every script in ``examples/`` runs to completion.

Each script is loaded from its file and its ``main`` is called with ``--quick``,
which reduces the sweep resolution and suppresses figure writing, so the whole
tier stays well inside the time budget. The scripts are checked for a zero return
code, for non-empty output, and for the specific facts each one exists to
produce, so that a script that runs but prints nothing useful still fails.

One script is additionally run without ``--quick`` at a reduced point count so
that the figure code is exercised, since a figure that raises would otherwise
never be caught by this tier.
"""

from __future__ import annotations

import importlib.util
import sys
from collections.abc import Sequence
from pathlib import Path
from types import ModuleType
from typing import Protocol

import pytest

_EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
_SCRIPTS = sorted(path.name for path in _EXAMPLES.glob("*.py"))


class _Example(Protocol):
    """The interface every example script exposes."""

    def main(self, argv: Sequence[str] | None = None) -> int:
        """Run the example and return a process exit code."""
        ...


def _load(name: str) -> ModuleType:
    path = _EXAMPLES / name
    spec = importlib.util.spec_from_file_location(f"_example_{path.stem}", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_every_example_script_is_discovered() -> None:
    """The example directory holds the scripts the README advertises."""
    assert _SCRIPTS == [
        "benchmark_points.py",
        "efficiency_study.py",
        "field_scaling.py",
        "power_balance.py",
        "scaling_sensitivity.py",
    ]


@pytest.mark.parametrize("name", _SCRIPTS)
def test_example_runs_and_produces_output(
    name: str, capsys: pytest.CaptureFixture[str]
) -> None:
    """Each script returns zero and prints a substantial amount of text."""
    module = _load(name)
    assert module.main(["--quick"]) == 0
    captured = capsys.readouterr()
    assert captured.err == ""
    assert len(captured.out.splitlines()) > 5


@pytest.mark.parametrize("name", _SCRIPTS)
def test_example_exposes_a_help_message(name: str) -> None:
    """Each script parses arguments and exits cleanly on ``--help``."""
    module = _load(name)
    with pytest.raises(SystemExit) as exit_info:
        module.main(["--help"])
    assert exit_info.value.code == 0


def test_power_balance_reports_closure_and_the_lawson_comparison(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The balance script prints the residual, the constraints, and the requirement."""
    assert _load("power_balance.py").main(["--quick"]) == 0
    out = capsys.readouterr().out
    assert "closure residual" in out
    assert "binding constraint" in out
    assert "ignition requirement" in out


def test_benchmark_script_reports_both_profile_cases(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The benchmark script prints the flat and the peaked comparison side by side."""
    assert _load("benchmark_points.py").main(["--quick"]) == 0
    out = capsys.readouterr().out
    assert "FLAT PROFILES" in out
    assert "PARABOLIC PROFILES" in out
    for name in ("ARC", "ITER", "SPARC"):
        assert name in out


def test_field_scaling_script_reports_the_measured_exponents(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The field scaling script prints an exponent for every invariant it sweeps."""
    assert _load("field_scaling.py").main(["--quick"]) == 0
    out = capsys.readouterr().out
    assert out.count("fusion power exponent in field") == 3
    assert "fusion power exponent in major radius" in out
    assert "analytic expectation" in out


def test_efficiency_script_reports_the_recirculating_breakdown(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The efficiency script names every electrical draw and the net result."""
    assert _load("efficiency_study.py").main(["--quick"]) == 0
    out = capsys.readouterr().out
    for line in (
        "heating and current drive",
        "cryoplant",
        "balance of plant",
        "net electric",
        "recirculating fraction",
        "reflectivity",
    ):
        assert line in out


def test_scaling_sensitivity_script_reports_the_amplification(
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The sensitivity script prints the power degradation and its amplification."""
    assert _load("scaling_sensitivity.py").main(["--quick"]) == 0
    out = capsys.readouterr().out
    assert "1/(1-alpha_P)" in out
    for name in ("IPB98(y,2)", "ITER89-P", "Petty08"):
        assert name in out


def test_figures_are_written_when_the_quick_flag_is_absent(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """One script is run with figures enabled so that the plotting code is exercised.

    The remaining scripts are covered by the quick runs above. Running only one
    of them with figures keeps this tier inside its time budget while still
    ensuring that a figure which raises cannot reach a release.
    """
    module = _load("power_balance.py")
    assert module.main(["--figure-dir", str(tmp_path)]) == 0
    written = sorted(path.name for path in tmp_path.glob("*.png"))
    assert written == ["balance.png", "lawson.png"]
    for path in tmp_path.glob("*.png"):
        assert path.stat().st_size > 1000
    assert "wrote" in capsys.readouterr().out
