#!/usr/bin/env bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <build directory> <output dir>"
    exit 1
fi

# Get arguments and resolve to absolute paths
buildDir="$(realpath "$1")"
outputDir="$(realpath "$2")"

# Get the directory where the script resides
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure output directory exists
mkdir -p "$outputDir"

# Create symlink to python script (absolute)
ln -sf "$script_dir/pyscript/test_semsolver.py" "$outputDir"

# Create symlinks to each python lib matching *cpython*
for lib in "$buildDir"/lib/*cpython*; do
    ln -sf "$lib" "$outputDir/"
done

# Create symlink to pykokkos base (absolute)
ln -sf "$buildDir/external/pykokkos-base/kokkos" "$outputDir/"
