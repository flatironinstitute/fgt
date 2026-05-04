"""High-level NumPy-friendly wrapper for the pfgt Fortran driver."""
from __future__ import annotations

import ctypes
from typing import Optional

import numpy as np

from . import _libfgt
from ._array import as_f64_fortran, as_int32, f64_ptr, zeros_f64


def pfgt(
    sources: np.ndarray,
    *,
    charges: Optional[np.ndarray] = None,
    rnormal: Optional[np.ndarray] = None,
    dipstr: Optional[np.ndarray] = None,
    targets: Optional[np.ndarray] = None,
    delta: float,
    eps: float = 1e-6,
    ifpgh: int = 1,
    ifpghtarg: int = 0,
    iperiod: int = 0,
    bs0: float = 1.0,
    cen0: Optional[np.ndarray] = None,
) -> dict:
    """Point fast Gauss transform.

    Parameters
    ----------
    sources : (dim, ns) float64
        Source locations. ``dim`` must be 1, 2, or 3.
    charges : (nd, ns) float64, optional
        Charge strengths. If None, no charge interactions.
    rnormal : (dim, ns) float64, optional
        Dipole directions. Required if ``dipstr`` is given.
    dipstr : (nd, ns) float64, optional
        Dipole strengths. If None, no dipole interactions.
    targets : (dim, nt) float64, optional
        Extra evaluation targets.
    delta : float
        Gaussian variance.
    eps : float, optional
        Requested precision (default 1e-6).
    ifpgh : int, optional
        At sources: 1 = pot, 2 = pot+grad, 3 = pot+grad+hess.
    ifpghtarg : int, optional
        At targets: 0 = none, 1, 2, or 3 (same meaning).
    iperiod : int, optional
        0 = free space (default), 1 = periodic.
    bs0 : float, optional
        Bounding-box size (default 1.0).
    cen0 : (dim,) float64, optional
        Bounding-box center (default = 0.5 * ones(dim)).

    Returns
    -------
    dict with keys 'pot', 'grad', 'hess', 'pottarg', 'gradtarg', 'hesstarg'
    where each entry is None if not requested.
    """
    sources = np.asarray(sources, dtype=np.float64)
    if sources.ndim != 2:
        raise ValueError(f"sources must be 2-D, got shape {sources.shape}")
    dim, ns = sources.shape
    if dim not in (1, 2, 3):
        raise ValueError(f"dim must be 1, 2, or 3, got {dim}")

    # nd = number of densities
    nd = 1
    if charges is not None:
        charges = np.asarray(charges, dtype=np.float64)
        if charges.ndim == 1:
            charges = charges[None, :]
        if charges.shape[1] != ns:
            raise ValueError(
                f"charges shape {charges.shape} mismatches ns={ns}")
        nd = charges.shape[0]

    if dipstr is not None:
        dipstr = np.asarray(dipstr, dtype=np.float64)
        if dipstr.ndim == 1:
            dipstr = dipstr[None, :]
        if dipstr.shape != (nd, ns):
            raise ValueError(
                f"dipstr shape {dipstr.shape} mismatches (nd={nd}, ns={ns})"
            )
        if rnormal is None:
            raise ValueError("rnormal must be provided when dipstr is given")
        rnormal = np.asarray(rnormal, dtype=np.float64)
        if rnormal.shape != (dim, ns):
            raise ValueError(
                f"rnormal shape {rnormal.shape} mismatches (dim={dim}, ns={ns})"
            )

    if cen0 is None:
        cen0 = np.full(dim, 0.5, dtype=np.float64)
    else:
        cen0 = np.asarray(cen0, dtype=np.float64)
        if cen0.shape != (dim,):
            raise ValueError(f"cen0 shape {cen0.shape} != (dim={dim},)")

    if targets is not None:
        targets = np.asarray(targets, dtype=np.float64)
        if targets.ndim != 2 or targets.shape[0] != dim:
            raise ValueError(
                f"targets must have shape (dim={dim}, nt), got {targets.shape}"
            )
        nt = targets.shape[1]
    else:
        nt = 0

    nhess = dim * (dim + 1) // 2

    # Allocate outputs in Fortran order.
    pot = zeros_f64((nd, max(ns, 1)))
    grad = zeros_f64((nd, dim, max(ns, 1)))
    hess = zeros_f64((nd, nhess, max(ns, 1)))
    pottarg = zeros_f64((nd, max(nt, 1)))
    gradtarg = zeros_f64((nd, dim, max(nt, 1)))
    hesstarg = zeros_f64((nd, nhess, max(nt, 1)))

    # Dummy filler arrays the Fortran side requires when the corresponding
    # flag is 0 (it still binds the dummy arg, even if it doesn't read it).
    arr_charge = charges if charges is not None else zeros_f64((nd, max(ns, 1)))
    arr_dipstr = dipstr if dipstr is not None else zeros_f64((nd, max(ns, 1)))
    arr_rnormal = rnormal if rnormal is not None else zeros_f64((dim, max(ns, 1)))
    arr_targ = targets if targets is not None else zeros_f64((dim, max(nt, 1)))

    cen0_f = as_f64_fortran(cen0, shape=(dim,), name="cen0")
    sources_f = as_f64_fortran(sources, shape=(dim, ns), name="sources")
    arr_charge = as_f64_fortran(arr_charge, name="charge")
    arr_dipstr = as_f64_fortran(arr_dipstr, name="dipstr")
    arr_rnormal = as_f64_fortran(arr_rnormal, name="rnormal")
    arr_targ = as_f64_fortran(arr_targ, name="targ")

    ifcharge = 1 if charges is not None else 0
    ifdipole = 1 if dipstr is not None else 0

    _libfgt.lib().fgt_pfgt_d(
        as_int32(nd, "nd"),
        as_int32(dim, "dim"),
        ctypes.c_double(delta),
        ctypes.c_double(eps),
        as_int32(iperiod, "iperiod"),
        ctypes.c_double(bs0),
        f64_ptr(cen0_f),
        as_int32(ns, "ns"),
        f64_ptr(sources_f),
        as_int32(ifcharge, "ifcharge"),
        f64_ptr(arr_charge),
        as_int32(ifdipole, "ifdipole"),
        f64_ptr(arr_rnormal),
        f64_ptr(arr_dipstr),
        as_int32(ifpgh, "ifpgh"),
        f64_ptr(pot),
        f64_ptr(grad),
        f64_ptr(hess),
        as_int32(nt, "nt"),
        f64_ptr(arr_targ),
        as_int32(ifpghtarg, "ifpghtarg"),
        f64_ptr(pottarg),
        f64_ptr(gradtarg),
        f64_ptr(hesstarg),
    )

    out: dict = {}
    out["pot"] = pot if ifpgh >= 1 else None
    out["grad"] = grad if ifpgh >= 2 else None
    out["hess"] = hess if ifpgh >= 3 else None
    out["pottarg"] = pottarg if ifpghtarg >= 1 else None
    out["gradtarg"] = gradtarg if ifpghtarg >= 2 else None
    out["hesstarg"] = hesstarg if ifpghtarg >= 3 else None
    return out
