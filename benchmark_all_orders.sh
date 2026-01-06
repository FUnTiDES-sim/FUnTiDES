#!/bin/bash
################################################################################
# Benchmark Script: Makutu vs Tensorial Implementation
# Tests all polynomial orders (Q1-Q5) with equivalent problem sizes
################################################################################

# Force C locale for consistent number formatting
export LC_ALL=C
export LANG=C

# Configuration
SEMPROXY="./build/bin/semproxy"
DOMAIN_SIZE="--lx 2000 --ly 2000 --lz 2000"
TIME_PARAMS="-timemax 1 --dt 0.001"
OUTPUT_DIR="benchmark_results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Summary file
SUMMARY_FILE="$OUTPUT_DIR/benchmark_summary_$TIMESTAMP.txt"

echo "================================================================================"
echo "  BENCHMARK: Makutu vs Tensorial Implementation (All Orders)"
echo "================================================================================"
echo "Timestamp: $(date)"
echo "Domain: 2000x2000x2000"
echo "Time: 0-1s with dt=0.001"
echo "Output: $OUTPUT_DIR"
echo "================================================================================"
echo ""

# Initialize summary
cat > "$SUMMARY_FILE" << EOF
BENCHMARK RESULTS - Makutu vs Tensorial
========================================
Date: $(date)
Domain: 2000x2000x2000
Time: 0-1s, dt=0.001

Order | Elements    | Makutu (s) | Tensorial (s) | Ratio   | Status
------|-------------|------------|---------------|---------|--------
EOF

################################################################################
# Function to run a single benchmark and extract time
################################################################################
run_benchmark() {
    local ORDER=$1
    local EX=$2
    local EY=$3
    local EZ=$4
    local IMPLEM=$5

    echo "Running: Order $ORDER, Grid ${EX}x${EY}x${EZ}, Implementation: $IMPLEM"

    OUTPUT_FILE="$OUTPUT_DIR/order${ORDER}_${IMPLEM}_$TIMESTAMP.log"

    # Run the benchmark
    $SEMPROXY --order $ORDER --ex $EX --ey $EY --ez $EZ \
              $DOMAIN_SIZE $TIME_PARAMS --implem $IMPLEM \
              > "$OUTPUT_FILE" 2>&1

    # Extract kernel time - format is "---- Elapsed Kernel Time : X.XXX seconds."
    KERNEL_TIME=$(grep "Elapsed Kernel Time" "$OUTPUT_FILE" | sed 's/.*: \(.*\) seconds.*/\1/')

    if [ -z "$KERNEL_TIME" ]; then
        echo "  ERROR: Could not extract kernel time"
        cat "$OUTPUT_FILE" | tail -20
        return 1
    fi

    echo "  -> Kernel Time: ${KERNEL_TIME}s"
    echo "$KERNEL_TIME"
}

################################################################################
# Function to compare results
################################################################################
compare_results() {
    local ORDER=$1
    local EX=$2
    local EY=$3
    local EZ=$4
    local MAKUTU_TIME=$5
    local TENSORIAL_TIME=$6

    # Calculate elements
    ELEMENTS=$((EX * EY * EZ))

    # Calculate ratio using awk
    RATIO=$(awk -v t="$TENSORIAL_TIME" -v m="$MAKUTU_TIME" 'BEGIN {printf "%.2f", t / m}')

    # Determine status
    STATUS="SLOWER"
    if awk -v r="$RATIO" 'BEGIN {exit !(r < 1.0)}'; then
        STATUS="FASTER"
    elif awk -v r="$RATIO" 'BEGIN {exit !(r < 1.1)}'; then
        STATUS="SIMILAR"
    fi

    # Print comparison
    printf "Order %d: Makutu=%.3fs, Tensorial=%.3fs, Ratio=%.2fx [%s]\n" \
           "$ORDER" "$MAKUTU_TIME" "$TENSORIAL_TIME" "$RATIO" "$STATUS"

    # Append to summary file
    printf "Q%-4d | %11d | %10.3f | %13.3f | %7.2fx | %s\n" \
           "$ORDER" "$ELEMENTS" "$MAKUTU_TIME" "$TENSORIAL_TIME" "$RATIO" "$STATUS" >> "$SUMMARY_FILE"
}

################################################################################
# Run Benchmarks for All Orders
################################################################################

echo ""
echo "================================================================================"
echo "  ORDER 1 (Q1) - Linear Elements"
echo "================================================================================"

tmp=$(run_benchmark 1 200 200 200 makutu )
MAKUTU_Q1=$(echo $tmp | awk '{print $12}')
tmp=$(run_benchmark 1 200 200 200 tensorial )
TENSORIAL_Q1=$(echo $tmp | awk '{print $12}')
compare_results 1 200 200 200 "$MAKUTU_Q1" "$TENSORIAL_Q1"
echo ""

echo "================================================================================"
echo "  ORDER 2 (Q2) - Quadratic Elements"
echo "================================================================================"
tmp=$(run_benchmark 2 100 100 100 makutu )
MAKUTU_Q2=$(echo $tmp | awk '{print $12}')
tmp=$(run_benchmark 2 100 100 100 tensorial )
TENSORIAL_Q2=$(echo $tmp | awk '{print $12}')
compare_results 2 100 100 100 "$MAKUTU_Q2" "$TENSORIAL_Q2"
echo ""

echo "================================================================================"
echo "  ORDER 3 (Q3) - Cubic Elements"
echo "================================================================================"
MAKUTU_Q3=$(run_benchmark 3 60 60 60 makutu)
TENSORIAL_Q3=$(run_benchmark 3 60 60 60 tensorial)
tmp=$(run_benchmark 3 60 60 100 makutu )
MAKUTU_Q3=$(echo $tmp | awk '{print $12}')
tmp=$(run_benchmark 3 60 60 100 tensorial )
TENSORIAL_Q3=$(echo $tmp | awk '{print $12}')
compare_results 3 60 60 60 "$MAKUTU_Q3" "$TENSORIAL_Q3"
echo ""

echo "================================================================================"
echo "  ORDER 4 (Q4) - Quartic Elements"
echo "================================================================================"
MAKUTU_Q4=$(run_benchmark 4 40 40 40 makutu)
TENSORIAL_Q4=$(run_benchmark 4 40 40 40 tensorial)
tmp=$(run_benchmark 4 40 40 40 makutu )
MAKUTU_Q4=$(echo $tmp | awk '{print $12}')
tmp=$(run_benchmark 4 40 40 40 tensorial )
TENSORIAL_Q4=$(echo $tmp | awk '{print $12}')
compare_results 4 40 40 40 "$MAKUTU_Q4" "$TENSORIAL_Q4"
echo ""

echo "================================================================================"
echo "  ORDER 5 (Q5) - Quintic Elements"
echo "================================================================================"
tmp=$(run_benchmark 5 30 30 30 makutu )
MAKUTU_Q5=$(echo $tmp | awk '{print $12}')
tmp=$(run_benchmark 5 30 30 30 tensorial )
TENSORIAL_Q5=$(echo $tmp | awk '{print $12}')
compare_results 5 30 30 30 "$MAKUTU_Q5" "$TENSORIAL_Q5"
echo ""

################################################################################
# Final Summary
################################################################################

echo "================================================================================"
echo "  BENCHMARK COMPLETE"
echo "================================================================================"
echo ""
cat "$SUMMARY_FILE"
echo ""
echo "Detailed logs saved to: $OUTPUT_DIR/"
echo "Summary saved to: $SUMMARY_FILE"
echo ""

# Create a simple CSV for plotting
CSV_FILE="$OUTPUT_DIR/benchmark_data_$TIMESTAMP.csv"
RATIO_Q1=$(awk -v t="$TENSORIAL_Q1" -v m="$MAKUTU_Q1" 'BEGIN {printf "%.3f", t / m}')
RATIO_Q2=$(awk -v t="$TENSORIAL_Q2" -v m="$MAKUTU_Q2" 'BEGIN {printf "%.3f", t / m}')
RATIO_Q3=$(awk -v t="$TENSORIAL_Q3" -v m="$MAKUTU_Q3" 'BEGIN {printf "%.3f", t / m}')
RATIO_Q4=$(awk -v t="$TENSORIAL_Q4" -v m="$MAKUTU_Q4" 'BEGIN {printf "%.3f", t / m}')
RATIO_Q5=$(awk -v t="$TENSORIAL_Q5" -v m="$MAKUTU_Q5" 'BEGIN {printf "%.3f", t / m}')

cat > "$CSV_FILE" << EOF
Order,Elements,Makutu_Time,Tensorial_Time,Ratio
1,8000000,$MAKUTU_Q1,$TENSORIAL_Q1,$RATIO_Q1
2,1000000,$MAKUTU_Q2,$TENSORIAL_Q2,$RATIO_Q2
3,216000,$MAKUTU_Q3,$TENSORIAL_Q3,$RATIO_Q3
4,64000,$MAKUTU_Q4,$TENSORIAL_Q4,$RATIO_Q4
5,27000,$MAKUTU_Q5,$TENSORIAL_Q5,$RATIO_Q5
EOF

echo "CSV data saved to: $CSV_FILE"
echo ""
echo "================================================================================"
