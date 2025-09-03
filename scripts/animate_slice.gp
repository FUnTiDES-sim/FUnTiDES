#!/bin/bash

echo "Starting automatic animation..."
echo "Press Ctrl+C to stop at any time"

# Create gnuplot script with automatic timing
cat > auto_animate.gp << 'EOF'
set terminal qt enhanced font "Arial,12"
set xlabel "X"
set ylabel "Y"
set size square
set palette defined (0 "blue", 0.5 "white", 1 "red")

# Get file count
n_files = system("ls slice*.dat | wc -l") + 0
print sprintf("Found %d slice files", n_files)

# Animation loop
do for [i=1:n_files] {
    filename = system(sprintf("ls slice*.dat | sort -V | sed -n '%dp'", i))

    # Read dimensions
    sizex = system(sprintf("head -1 %s", filename)) + 0
    sizey = system(sprintf("head -2 %s | tail -1", filename)) + 0

    print sprintf("Frame %d/%d: %s", i, n_files, filename)

    # Set ranges
    set xrange [0:sizex-1]
    set yrange [0:sizey-1]
    set title sprintf("Frame %d/%d: %s (%dx%d)", i, n_files, filename, sizex, sizey)

    # Get data range and center colorbar on 0
    stats filename skip 2 matrix nooutput
    max_abs = (abs(STATS_max) > abs(STATS_min)) ? abs(STATS_max) : abs(STATS_min)
    set cbrange [-max_abs:max_abs]

    # Plot with centered colorbar
    plot filename skip 2 matrix with image notitle

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
