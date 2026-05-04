# fgt — Python interface

Point and box fast Gauss transforms in 1, 2, 3 dimensions.

This is a thin `ctypes` wrapper over the Fortran library `libfgt`.

```python
import numpy as np
import fgt

rng = np.random.default_rng(0)
sources = rng.uniform(size=(3, 1000))
charges = rng.standard_normal((1, 1000))
out = fgt.pfgt(sources, charges=charges, delta=1e-3, eps=1e-6, ifpgh=1)
```

See the top-level project [README](https://github.com/sj90101/fgt) for full
documentation.
