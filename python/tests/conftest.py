"""pytest fixtures for fgt tests."""
from __future__ import annotations

import pytest


@pytest.fixture(autouse=True)
def _ensure_lib_loadable():
    """Ensure libfgt is loadable; skip the test if not."""
    try:
        from fgt import _libfgt  # noqa: F401
    except OSError as e:
        pytest.skip(f"libfgt not loadable: {e}")
