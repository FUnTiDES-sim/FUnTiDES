#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import glob
import re
import os
import subprocess
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import MaxNLocator

# --- Configuration ---
frame_dir = "frames"
framerate = 10
video_filename = "heatmap_video.mp4"
colormap = "PuOr"
linthresh = 1e-2  # fixed linear threshold around zero
dpi = 150         # output image resolution
figsize = (10, 8) # size of each frame in inches

# --- Setup ---
os.makedirs(frame_dir, exist_ok=True)

def extract_number(filename):
    match = re.search(r"slice(\d+)\.dat", filename)
    return int(match.group(1)) if match else -1

# --- Step 1: Find all files and compute global min/max for colormap scaling ---
files = sorted(glob.glob("slice*.dat"), key=extract_number)

vmin, vmax = float("inf"), float("-inf")
matrix_shapes = None

for file in files:
    with open(file) as f:
        lines = f.readlines()
        rows = int(lines[0])
        cols = int(lines[1])
        data = list(map(float, lines[2].split()))
        matrix = np.array(data).reshape((rows, cols))

        if matrix_shapes is None:
            matrix_shapes = (rows, cols)
        elif matrix.shape != matrix_shapes:
            raise ValueError(f"Inconsistent matrix shape in {file}")

        vmin = min(vmin, matrix.min())
        vmax = max(vmax, matrix.max())

print(f"Global min: {vmin}, max: {vmax}, linthresh: {linthresh}")

# --- Prepare normalization of colormap ---
norm = SymLogNorm(linthresh=linthresh, vmin=vmin, vmax=vmax, base=10)

# --- Step 2: Generate annotated heatmap frames ---
for idx, file in enumerate(files):
    with open(file) as f:
        lines = f.readlines()
        rows = int(lines[0])
        cols = int(lines[1])
        data = list(map(float, lines[2].split()))
        matrix = np.array(data).reshape((rows, cols))

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(matrix, cmap=colormap, norm=norm, interpolation='nearest', aspect='auto')

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))

    plt.colorbar(im, ax=ax)
    ax.set_title(f"Frame {idx}")
    plt.savefig(f"{frame_dir}/frame_{idx:04d}.png", dpi=dpi)
    plt.close()

# --- Step 3: Create video with ffmpeg ---
subprocess.run([
    "ffmpeg", "-y", "-framerate", str(framerate),
    "-i", f"{frame_dir}/frame_%04d.png",
    "-c:v", "libx264",
    "-preset", "slow",
    "-crf", "18",
    "-pix_fmt", "yuv420p",
    video_filename
])

print(f"Video saved as {video_filename}")
