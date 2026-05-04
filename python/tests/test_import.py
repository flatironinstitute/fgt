"""Confirm fgt loads and exposes the expected symbols."""
from __future__ import annotations


def test_import_version():
    import fgt

    assert isinstance(fgt.__version__, str)
    assert "." in fgt.__version__


def test_pfgt_callable():
    import fgt

    assert callable(fgt.pfgt)


def test_libfgt_loaded():
    from fgt import _libfgt

    lib = _libfgt.lib()
    assert hasattr(lib, "fgt_pfgt_d")
    assert hasattr(lib, "fgt_bfgt_d")
