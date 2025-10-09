#!/bin/bash

# Script to compare benchmark files between baseline and latest results

if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <baseline_dir> <comparison_dir>"
  exit 1
fi

BASELINE_DIR="$1"
COMPARISON_DIR="$2"

if [ ! -d "$BASELINE_DIR" ]; then
  echo "Error: Baseline directory '$BASELINE_DIR' does not exist"
  exit 1
fi

if [ ! -d "$COMPARISON_DIR" ]; then
  echo "Error: Comparison directory '$COMPARISON_DIR' does not exist"
  exit 1
fi

# Find all JSON files in baseline directory
find "$BASELINE_DIR" -type f -name "*.json" | while read -r baseline_file; do
  # Get relative path from baseline directory
  rel_path="${baseline_file#$BASELINE_DIR/}"
  
  # Get just the filename
  filename=$(basename "$rel_path")
  dirname=$(dirname "$rel_path")
  
  # Transform baseline filename to comparison filename format
  # Remove extension, add bench_ prefix and _latest suffix
  base_name="${filename%.json}"
  comparison_filename="bench_${base_name}_latest.json"
  
  # Construct corresponding comparison file path
  if [ "$dirname" = "." ]; then
    comparison_file="$COMPARISON_DIR/$comparison_filename"
  else
    comparison_file="$COMPARISON_DIR/$dirname/$comparison_filename"
  fi
  
  if [ -f "$comparison_file" ]; then
    echo "Comparing: $filename -> $comparison_filename"
    python ./external/benchmark/tools/compare.py benchmarks "$baseline_file" "$comparison_file"
    echo "---"
  else
    echo "Warning: No matching file found for $filename (expected: $comparison_filename)"
  fi
done