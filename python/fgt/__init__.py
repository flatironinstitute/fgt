"""fgt - point and box fast Gauss transforms in 1, 2, 3 dimensions."""
from __future__ import annotations

__version__ = "0.1.0"

from .pfgt import pfgt

__all__ = ["pfgt", "__version__"]
