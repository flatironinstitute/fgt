#!/usr/bin/env bash
# python/scripts/build-libfgt.sh
#
# Build libfgt + bundled finufft from a fresh checkout. Invoked by
# cibuildwheel before-build on each platform.
#
# Honors environment overrides FC, CC, CXX, FFTW_ROOT, FGT_MAKEINC.

set -euo pipefail

cd "$(dirname "$0")/../.."

echo "[build-libfgt.sh] CWD=$(pwd)"
echo "[build-libfgt.sh] uname: $(uname -s) $(uname -m)"

# On manylinux_2_28 we install gcc-toolset-13 in before-all; activate it
# now so all gcc/gfortran/ld invocations use the toolset binaries.
if [ -f /opt/rh/gcc-toolset-13/enable ]; then
    set +u
    # shellcheck disable=SC1091
    source /opt/rh/gcc-toolset-13/enable
    set -u
    echo "[build-libfgt.sh] activated gcc-toolset-13: $(gcc --version | head -1)"
fi

# Initialize submodules if the runner forgot --recurse-submodules.
git submodule update --init --recursive

# Pick a make.inc appropriate for the platform (overridable via FGT_MAKEINC).
case "$(uname -s)" in
    Darwin) MAKEINC=${FGT_MAKEINC:-make.inc.macos.gnu} ;;
    Linux)  MAKEINC=${FGT_MAKEINC:-make.inc.linux.gnu.openblas} ;;
    MINGW*|MSYS*|CYGWIN*) MAKEINC=${FGT_MAKEINC:-make.inc.windows.mingw} ;;
    *) MAKEINC=${FGT_MAKEINC:-make.inc} ;;
esac
echo "[build-libfgt.sh] using $MAKEINC"
cp "$MAKEINC" make.inc

# Always clean: cibuildwheel mounts the host filesystem, so stale .o files
# from a prior host build (potentially with a different ABI) would otherwise
# be reused. Same logic for the finufft submodule.
make objclean || true
rm -f lib/*.so lib/*.dylib lib/*.dll lib-static/*.a 2>/dev/null || true
make -C external/finufft objclean || true
rm -f external/finufft/lib/*.so external/finufft/lib/*.dylib \
      external/finufft/lib-static/*.a 2>/dev/null || true

# Build (uses opts%nthreads=1 from fgt source; finufft is OMP=OFF + LTO=OFF
# via the fgt makefile's finufft target). LTO is disabled in the finufft
# submake itself (LTO=OFF), so passing -fno-lto on the make command line is
# unnecessary and previously clobbered finufft's own CXXFLAGS, which carry
# the -Iinclude -Ideps/xsimd/include paths.
make lib OMP=ON -j 4

ls -lh lib lib-static
