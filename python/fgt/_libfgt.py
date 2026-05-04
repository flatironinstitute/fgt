"""Locate and load libfgt via ctypes, and declare function prototypes.

The shared library ships either inside this package (``fgt/lib/``) for
wheel installs, or is found on the system loader path for from-source
installs (``DYLD_LIBRARY_PATH`` / ``LD_LIBRARY_PATH``).
"""
from __future__ import annotations

import ctypes
import os
import sys
from importlib import resources
from pathlib import Path

_LIBNAME_CANDIDATES = {
    "linux": ["libfgt.so"],
    "darwin": ["libfgt.so", "libfgt.dylib"],
    "win32": ["libfgt.dll", "fgt.dll"],
}


def _platform_key() -> str:
    if sys.platform.startswith("linux"):
        return "linux"
    if sys.platform == "darwin":
        return "darwin"
    if sys.platform == "win32":
        return "win32"
    raise RuntimeError(f"Unsupported platform: {sys.platform}")


def _candidate_paths() -> list[Path]:
    """Where to look for libfgt, in priority order."""
    out: list[Path] = []
    # 1) bundled with the wheel: <package>/lib/
    try:
        with resources.as_file(resources.files("fgt").joinpath("lib")) as p:
            out.append(Path(p))
    except (FileNotFoundError, ModuleNotFoundError):
        pass
    # 2) repo-local dev install: ../../lib relative to this file
    out.append(Path(__file__).resolve().parents[2] / "lib")
    # 3) typical install prefixes
    if "FGT_LIB_DIR" in os.environ:
        out.append(Path(os.environ["FGT_LIB_DIR"]))
    out.append(Path.home() / "lib")
    out.append(Path("/usr/local/lib"))
    out.append(Path("/usr/lib"))
    return out


def _load_libfgt() -> ctypes.CDLL:
    names = _LIBNAME_CANDIDATES[_platform_key()]
    last_err: Exception | None = None
    for d in _candidate_paths():
        for n in names:
            p = d / n
            if p.exists():
                try:
                    return ctypes.CDLL(str(p))
                except OSError as e:
                    last_err = e
                    continue
    msg = (
        f"Could not locate libfgt on this system. Searched: "
        f"{[str(p) for p in _candidate_paths()]}"
    )
    if last_err is not None:
        msg += f"\nLast loader error: {last_err}"
    raise OSError(msg)


# ----------------------------------------------------------------------
# Prototypes
# ----------------------------------------------------------------------
_lib = _load_libfgt()

_DBL_P = ctypes.POINTER(ctypes.c_double)
_INT_P = ctypes.POINTER(ctypes.c_int)


_lib.fgt_pfgt_d.restype = None
_lib.fgt_pfgt_d.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
    ctypes.c_int, ctypes.c_double, _DBL_P,
    ctypes.c_int, _DBL_P,
    ctypes.c_int, _DBL_P, ctypes.c_int, _DBL_P, _DBL_P,
    ctypes.c_int, _DBL_P, _DBL_P, _DBL_P,
    ctypes.c_int, _DBL_P,
    ctypes.c_int, _DBL_P, _DBL_P, _DBL_P,
]


_lib.fgt_bfgt_d.restype = None
_lib.fgt_bfgt_d.argtypes = [
    ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double,
    ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_int,
    ctypes.c_int, ctypes.c_int, ctypes.c_int,
    _INT_P, _INT_P,
    _DBL_P, _DBL_P, _DBL_P,
    ctypes.c_int, _DBL_P, _DBL_P, _DBL_P,
    ctypes.c_int,
    ctypes.c_int, _DBL_P, ctypes.c_int,
    _DBL_P, _DBL_P, _DBL_P,
    _DBL_P,
]


def lib() -> ctypes.CDLL:
    """Return the loaded libfgt CDLL handle (for tests)."""
    return _lib
