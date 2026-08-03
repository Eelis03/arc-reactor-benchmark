"""Tier one: what an installed copy of this package has to carry.

Passing mypy in strict mode inside this repository says nothing about what the
package delivers to anything that installs it. Without the PEP 561 marker file, a
downstream type checker treats every symbol here as untyped and silently skips
it, and no other test in this suite would notice. That is what this file exists
to catch.
"""

from __future__ import annotations

from pathlib import Path

import arc_benchmark

_MARKER = "py.typed"


def _package_directory() -> Path:
    """The directory the package is imported from, wherever it is installed."""
    assert arc_benchmark.__file__ is not None, "arc_benchmark must be a regular package"
    return Path(arc_benchmark.__file__).resolve().parent


def test_the_typing_marker_sits_inside_the_package_directory() -> None:
    """The marker is beside ``__init__.py``, which is where PEP 561 requires it.

    Checked against the imported location rather than against a path built from
    this file, so that an installed copy missing the marker fails here even
    though the marker is present in the source tree.
    """
    package = _package_directory()
    assert package.name == "arc_benchmark"
    assert (package / "__init__.py").is_file()
    assert (package / _MARKER).is_file()


def test_the_typing_marker_is_empty() -> None:
    """PEP 561 reads only the presence of the file, so it carries nothing.

    An empty file makes that explicit. Content here would be content nothing
    parses.
    """
    assert (_package_directory() / _MARKER).read_bytes() == b""
