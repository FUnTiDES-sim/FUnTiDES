# Benchmark Configuration

## System overview
 
Benchmarks were generated on Plafrim's sirocco09 machine:
 
- Processor: 2x 16 Core Intel Broadwell 
- Memory: 256GB per node
- GPU: 2 Nvidia P100

## Compiler and runtime environment

Loaded Modulefiles:
  1) tools/git/2.36.0      2) compiler/cuda/12.3    3) build/cmake/3.26.0    4) compiler/gcc/11.2.0

## Installation

Installed with Funtides-TPL.

Commit used is:
22601e40f81c07c86a3fc8fa6fef43416fafb01d


## Benchmark execution

Using benchmark bash script for GPU/Laptop. Node was reserved using slurm and `--exclusive` option.
Benchmarks were launched in interactive mode.
