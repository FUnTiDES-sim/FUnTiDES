#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

source "$PROJECT_ROOT/scripts/environments/env_Maple_GH200.sh"

BUILD_DIR="$PROJECT_ROOT/build"

# Initialize results file
RESULTS_FILE="$BUILD_DIR/timing_results.md"
> "$RESULTS_FILE"

# Write markdown table header
echo "| Type     | Grid       | X   | Y   | Z   | Param | Initial  | Compute | Total   |" >> "$RESULTS_FILE"
echo "|----------|------------|-----|-----|-----|-------|----------|---------|---------|" >> "$RESULTS_FILE"

# Function to run benchmark and extract timing
run_benchmark() {
    local solver=$1
    local mesh=$2
    local ex=$3
    local ey=$4
    local ez=$5
    local order=$6
    
    echo "Running: solver=$solver, mesh=$mesh, ex=$ex, ey=$ey, ez=$ez, order=$order"
    
    if [ "$solver" = "elastic" ]; then
        OUTPUT=$($BUILD_DIR/bin/semproxy --ex $ex --ey $ey --ez $ez --method=sem --implem=makutu --mesh=$mesh -o $order --dt 0.000963391 --timemax 1.5 --is-elastic 2>&1)
    else
        OUTPUT=$($BUILD_DIR/bin/semproxy --ex $ex --ey $ey --ez $ez --method=sem --implem=makutu --mesh=$mesh -o $order --dt 0.000963391 --timemax 1.5 2>&1)
    fi
    
    INIT_TIME=$(echo "$OUTPUT" | grep "Elapsed Initial Time" | awk '{print $5}')
    COMPUTE_TIME=$(echo "$OUTPUT" | grep "Elapsed Compute Time" | awk '{print $5}')
    TOTAL_TIME=$(echo "$OUTPUT" | grep "Elapsed TotalExe Time" | awk '{print $5}')
    
    echo "| $solver | $mesh | $ex | $ey | $ez | $order | $INIT_TIME | $COMPUTE_TIME | $TOTAL_TIME |" >> "$RESULTS_FILE"
}

# Elastic Solver - cartesian mesh
echo "=== Elastic Solver - Cartesian Mesh ==="
run_benchmark elastic cartesian 400 400 400 1
run_benchmark elastic cartesian 200 200 200 2
run_benchmark elastic cartesian 133 133 133 3

# Acoustic Solver - cartesian mesh
echo "=== Acoustic Solver - Cartesian Mesh ==="
run_benchmark acoustic cartesian 400 400 400 1
run_benchmark acoustic cartesian 200 200 200 2
run_benchmark acoustic cartesian 133 133 133 3

# Elastic Solver - ucartesian mesh
echo "=== Elastic Solver - Unstructured Cartesian Mesh ==="
run_benchmark elastic ucartesian 400 400 400 1
run_benchmark elastic ucartesian 200 200 200 2
run_benchmark elastic ucartesian 133 133 133 3

# Acoustic Solver - ucartesian mesh
echo "=== Acoustic Solver - Unstructured Cartesian Mesh ==="
run_benchmark acoustic ucartesian 400 400 400 1
run_benchmark acoustic ucartesian 200 200 200 2
run_benchmark acoustic ucartesian 133 133 133 3

# Display results in table format
echo ""
echo "==================== TIMING RESULTS ===================="
echo ""
printf "%-10s %-12s %-10s %-8s %-12s %-12s %-12s\n" "Solver" "Mesh" "Elements" "Order" "Init(s)" "Compute(s)" "Total(s)"
printf "%-10s %-12s %-10s %-8s %-12s %-12s %-12s\n" "----------" "------------" "----------" "--------" "------------" "------------" "------------"

while IFS=',' read -r solver mesh ex ey ez order init compute total; do
    elements="${ex}x${ey}x${ez}"
    printf "%-10s %-12s %-10s %-8s %-12s %-12s %-12s\n" "$solver" "$mesh" "$elements" "$order" "$init" "$compute" "$total"
done < "$RESULTS_FILE"

echo ""
echo "========================================================"
echo ""
echo "Detailed results saved to: $RESULTS_FILE"