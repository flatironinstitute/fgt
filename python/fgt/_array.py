"""NumPy <-> Fortran-contiguous helpers."""
from __future__ import annotations

import ctypes
from typing import Optional

import numpy as np


def as_f64_fortran(a: Optional[np.ndarray], shape: tuple[int, ...] | None = None,
                   name: str = "array") -> Optional[np.ndarray]:
    """Return a Fortran-contiguous float64 view of *a*, copying if needed.

    If *a* is None, return None.
    If *shape* is given, validates a.shape == shape.
    """
    if a is None:
        return None
    arr = np.ascontiguousarray(a, dtype=np.float64)
    arr = np.asfortranarray(arr)
    if shape is not None and arr.shape != shape:
        raise ValueError(f"{name}: expected shape {shape}, got {arr.shape}")
    return arr


def as_int32(x: int | np.ndarray, name: str = "scalar") -> int:
    """Return *x* as a Python int, validating it fits in int32."""
    v = int(x)
    if not (-(2**31) <= v < 2**31):
        raise OverflowError(f"{name}={v} does not fit in int32")
    return v


def f64_ptr(a: Optional[np.ndarray]):
    """Return a ctypes pointer to a's data, or NULL for None."""
    if a is None:
        return ctypes.cast(ctypes.c_void_p(0),
                           ctypes.POINTER(ctypes.c_double))
    if not (a.flags["F_CONTIGUOUS"] and a.dtype == np.float64):
        raise TypeError("f64_ptr requires float64 + Fortran-contiguous array")
    return a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def empty_f64(shape: tuple[int, ...]) -> np.ndarray:
    """Allocate a fresh Fortran-contiguous float64 array with given shape."""
    return np.empty(shape, dtype=np.float64, order="F")


def zeros_f64(shape: tuple[int, ...]) -> np.ndarray:
    return np.zeros(shape, dtype=np.float64, order="F")
