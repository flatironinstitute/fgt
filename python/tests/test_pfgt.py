"""Compare fgt.pfgt against direct O(N^2) Gauss summation."""
from __future__ import annotations

import numpy as np
import pytest

import fgt
from ._direct import pfgt_direct


@pytest.mark.parametrize("dim", [2, 3])
def test_pfgt_charges_only_pot(dim: int) -> None:
    rng = np.random.default_rng(42)
    ns = 200
    eps = 1e-6
    delta = 1e-2

    sources = np.asfortranarray(rng.uniform(0.1, 0.9, size=(dim, ns)))
    charges = np.asfortranarray(rng.standard_normal((1, ns)))

    out = fgt.pfgt(
        sources, charges=charges,
        delta=delta, eps=eps,
        ifpgh=1, ifpghtarg=0,
    )

    ref_src, _ = pfgt_direct(sources, charges, None, delta)
    err = np.max(np.abs(out["pot"][:, :ns] - ref_src))
    assert err < 2e-5, f"err={err}"


@pytest.mark.parametrize("dim", [2, 3])
def test_pfgt_targets(dim: int) -> None:
    rng = np.random.default_rng(0)
    ns, nt = 200, 50
    eps = 1e-6
    delta = 1e-2

    sources = np.asfortranarray(rng.uniform(0.1, 0.9, size=(dim, ns)))
    charges = np.asfortranarray(rng.standard_normal((1, ns)))
    targets = np.asfortranarray(rng.uniform(0.1, 0.9, size=(dim, nt)))

    out = fgt.pfgt(
        sources, charges=charges, targets=targets,
        delta=delta, eps=eps,
        ifpgh=1, ifpghtarg=1,
    )

    _, ref_tgt = pfgt_direct(sources, charges, targets, delta)
    err = np.max(np.abs(out["pottarg"][:, :nt] - ref_tgt))
    assert err < 2e-5, f"err={err}"
