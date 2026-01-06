#!/bin/bash
################################################################################
# Quick Benchmark: Makutu vs Tensorial (All Orders Q1-Q5)
################################################################################

echo "================================================================================"
echo "  TENSORIAL vs MAKUTU BENCHMARK"
echo "================================================================================"
echo "Start time: $(date)"
echo ""

# Order 1
echo "--- ORDER 1 (200x200x200) ---"
./semproxy --order 1 --ex 200 --ey 200 --ez 200 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem makutu 2>&1 | grep "Elapsed Kernel Time"
./semproxy --order 1 --ex 200 --ey 200 --ez 200 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem tensorial 2>&1 | grep "Elapsed Kernel Time"
echo ""

# Order 2
echo "--- ORDER 2 (100x100x100) ---"
./semproxy --order 2 --ex 100 --ey 100 --ez 100 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem makutu 2>&1 | grep "Elapsed Kernel Time"
./semproxy --order 2 --ex 100 --ey 100 --ez 100 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem tensorial 2>&1 | grep "Elapsed Kernel Time"
echo ""

# Order 3
echo "--- ORDER 3 (60x60x60) ---"
./semproxy --order 3 --ex 60 --ey 60 --ez 60 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem makutu 2>&1 | grep "Elapsed Kernel Time"
./semproxy --order 3 --ex 60 --ey 60 --ez 60 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem tensorial 2>&1 | grep "Elapsed Kernel Time"
echo ""

# Order 4
echo "--- ORDER 4 (40x40x40) ---"
./semproxy --order 4 --ex 40 --ey 40 --ez 40 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem makutu 2>&1 | grep "Elapsed Kernel Time"
./semproxy --order 4 --ex 40 --ey 40 --ez 40 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem tensorial 2>&1 | grep "Elapsed Kernel Time"
echo ""

# Order 5
echo "--- ORDER 5 (30x30x30) ---"
./semproxy --order 5 --ex 30 --ey 30 --ez 30 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem makutu 2>&1 | grep "Elapsed Kernel Time"
./semproxy --order 5 --ex 30 --ey 30 --ez 30 --lx 2000 --ly 2000 --lz 2000 -timemax 1 --dt 0.001 --implem tensorial 2>&1 | grep "Elapsed Kernel Time"
echo ""

echo "================================================================================"
echo "End time: $(date)"
echo "================================================================================"
