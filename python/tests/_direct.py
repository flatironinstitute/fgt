"""Direct O(N^2) evaluation of the Gauss transform - reference for tests.

For small N this is the ground truth we compare against.
"""
from __future__ import annotations

import numpy as np


def pfgt_direct(
    sources: np.ndarray,
    charges: np.ndarray | None,
    targets: np.ndarray | None,
    delta: float,
):
    """Compute pot at sources and (optionally) at targets via direct sum.

    pot[k, i] = sum_j charges[k, j] * exp(-||x_i - x_j||^2 / delta)

    Self-interaction (i == j) is included.
    """
    if charges is None:
        return None, None

    diff = sources[:, :, None] - sources[:, None, :]   # (dim, ns, ns)
    sq = np.einsum("dij,dij->ij", diff, diff)
    K = np.exp(-sq / delta)
    pot_src = charges @ K                               # (nd, ns)

    pot_tgt = None
    if targets is not None:
        diff_t = targets[:, :, None] - sources[:, None, :]   # (dim, nt, ns)
        sq_t = np.einsum("dij,dij->ij", diff_t, diff_t)
        K_t = np.exp(-sq_t / delta)
        pot_tgt = charges @ K_t.T                            # (nd, nt)

    return pot_src, pot_tgt
