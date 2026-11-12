# Environment to build and test a FUnTiDES build tree that shares its TPL install with an installed FUnTiDES.
#
# Usage (from any shell):
#   env -u LD_LIBRARY_PATH bash -c 'source scripts/dev/build_tree_env.sh build-2d; cmake --build build-2d -j8'
#
# Why: the TPL environment (<install>/setup_env.sh) prepends <install>/lib to LD_LIBRARY_PATH, and
# LD_LIBRARY_PATH is searched before the RUNPATH of the test executables. When FUnTiDES itself is installed
# in the same prefix, the tests then load the installed libfun_*.so instead of the build tree's: as soon
# as a header adds a symbol, gtest discovery fails at build time with "undefined symbol". This script
# sources the TPL environment, then puts the build tree's library directories first.
#
# $1: build directory (default: build). TPL environment: $FUNTIDES_TPL_ENV, else <Kokkos_DIR>/../../../setup_env.sh.

_bt_dir="$(cd "${1:-build}" && pwd)"
_bt_kokkos="$(sed -n 's/^Kokkos_DIR:[A-Z]*=//p' "${_bt_dir}/CMakeCache.txt" | head -1)"
_bt_tpl="${FUNTIDES_TPL_ENV:-${_bt_kokkos%/lib/cmake/Kokkos}/setup_env.sh}"
if [ -f "${_bt_tpl}" ]; then
    # shellcheck disable=SC1090
    source "${_bt_tpl}" >/dev/null 2>&1
else
    echo "build_tree_env.sh: TPL environment not found (${_bt_tpl}); set FUNTIDES_TPL_ENV" >&2
fi
_bt_libs="$(find "${_bt_dir}/src" "${_bt_dir}/lib" -name 'libfun_*.so' -printf '%h\n' 2>/dev/null | sort -u | paste -sd:)"
export LD_LIBRARY_PATH="${_bt_libs}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
unset _bt_dir _bt_kokkos _bt_tpl _bt_libs
