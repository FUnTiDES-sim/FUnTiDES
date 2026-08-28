#!/bin/bash

echo "Starting automatic animation..."
echo "Press Ctrl+C to stop at any time"

# Create gnuplot script with automatic timing
cat > auto_animate.gp << 'EOF'
set terminal qt enhanced font "Arial,12"
set xlabel "X"
set ylabel "Y"
set size square
set palette defined (0 "black", 0.5 "gray", 1 "white")

# Get file count
n_files = system("ls slice*.dat | wc -l") + 0
print sprintf("Found %d slice files", n_files)

# Animation loop
do for [i=1:n_files] {
    filename = system(sprintf("ls slice*.dat | sort -V | sed -n '%dp'", i))

    print sprintf("Frame %d/%d: %s", i, n_files, filename)

    # Get data range and dimensions in one pass (no header lines in these files)
    stats filename matrix nooutput

    set xrange [0:STATS_size_x-1]
    set yrange [0:STATS_size_y-1]
    set title sprintf("Frame %d/%d: %s (%dx%d)", i, n_files, filename, STATS_size_x, STATS_size_y)

    # Check if data has any variation
    if (STATS_max == STATS_min) {
        # All values are the same (e.g., all zeros)
        if (abs(STATS_max) > 0) {
            # All same non-zero value
            set cbrange [STATS_max*0.9:STATS_max*1.1]
        } else {
            # All zeros - use small range
            set cbrange [-1e-10:1e-10]
        }
        print sprintf("  Warning: All values are %.6f", STATS_max)
    } else {
        # Normal case - center colorbar on 0
        max_abs = (abs(STATS_max) > abs(STATS_min)) ? abs(STATS_max) : abs(STATS_min)
        set cbrange [-max_abs:max_abs]
    }

    # Plot with centered colorbar
    plot filename matrix with image notitle

    # Small delay between frames
    pause 0.15
}

print "Animation complete"
pause -1
EOF

# Run the animation
gnuplot auto_animate.gp

# Cleanup
rm -f auto_animate.gp
