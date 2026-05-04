"""setup.py - thin wrapper that copies a prebuilt libfgt.{so,dylib,dll}
into the wheel.

The actual library is built by ``python/scripts/build-libfgt.sh`` (invoked
by cibuildwheel's ``before-build``) before this script runs. We just
bundle whatever shared library lives at ./lib/.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution

try:
    from wheel.bdist_wheel import bdist_wheel
except ImportError:
    bdist_wheel = None


HERE = Path(__file__).resolve().parent  # repo root, where pyproject.toml lives
LIB_SRC = HERE / "lib"
LIB_DST = HERE / "python" / "fgt" / "lib"


def _libfgt_filename() -> str:
    if sys.platform.startswith("linux"):
        return "libfgt.so"
    if sys.platform == "darwin":
        for cand in ("libfgt.so", "libfgt.dylib"):
            if (LIB_SRC / cand).exists():
                return cand
        return "libfgt.so"
    if sys.platform == "win32":
        return "libfgt.dll"
    raise RuntimeError(f"Unsupported platform: {sys.platform}")


class build_py_with_lib(build_py):
    """Copy libfgt.{so,dylib,dll} into fgt/lib/ before packaging."""

    def run(self) -> None:
        LIB_DST.mkdir(parents=True, exist_ok=True)
        name = _libfgt_filename()
        src = LIB_SRC / name
        if not src.exists():
            raise RuntimeError(
                f"Prebuilt {src} not found. Run the build script first:\n"
                f"  bash {HERE}/python/scripts/build-libfgt.sh"
            )
        dst = LIB_DST / name
        shutil.copy2(src, dst)
        print(f"copied {src} -> {dst}")
        super().run()


class BinaryDistribution(Distribution):
    """Force the wheel to be marked platform-specific (it contains a .so)."""

    def has_ext_modules(self):
        return True


cmdclass = {"build_py": build_py_with_lib}

if bdist_wheel is not None:

    class bdist_wheel_py3_plat(bdist_wheel):
        """Tag as py3-none-<plat>: ctypes uses no cpython-API symbols."""

        def get_tag(self):
            _python, _abi, plat = super().get_tag()
            return ("py3", "none", plat)

    cmdclass["bdist_wheel"] = bdist_wheel_py3_plat


setup(
    cmdclass=cmdclass,
    distclass=BinaryDistribution,
)
