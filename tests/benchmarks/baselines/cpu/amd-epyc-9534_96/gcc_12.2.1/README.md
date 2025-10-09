# Benchmark Configuration

## System overview

CHECK WHAT IS PUBLIC


## Compiler and runtime environment

- Programming Environment: `PrgEnv-gnu/8.4.0`
- Cray Programming Environment: `craype/2.7.23`
- Compiler: `gcc/12.2.1`
- CPU Target Module: `craype-x86-milan`
- Python Version: `python/3.11.6`
- NumPy Version: `numpy/2.3.3`

## Build configuration

```bash
source scripts/environments/env_Pangea4.sh
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
pip install -r requirements-dev.txt
mkdir build
cd build
cmake -DUSE_KOKKOS=ON -DENABLE_PYWRAP=ON -DCMAKE_INSTALL_PREFIX=../install  ..
make install
```

## Benchmark execution
Benchmarks were executed using the SLURM batch script [run_benchmarks_p4.sbatch](../../../../../../scripts/benchmarks/run_benchmarks_p4.sbatch).
```bash
sbatch scripts/benchmarks/run_benchmarks_p4.sbatch
```

## Job configuration
Slurm settings:
- Partition: p4_dev
- Resources: 1 full socket (96 cores), exclusive node access
- Memory: 734GB
- Threading: Hyperthreading disabled, cores bound to socket

OpenMP settings:
- `OMP_PLACES=cores`
- `OMP_PROC_BIND=close`
- `OMP_DYNAMIC=FALSE`

Benchmark suites:
- Python Benchmarks: run_pywrap_benchmarks.py
- C++ Benchmarks: CTest suite (labeled benchmark)

## Generate plots

To plot current baselines
```bash
python scripts/benchmarks/plot_benchmarks.py --python-dir tests/benchmarks/baselines/cpu/amd-epyc-9534_96/gcc_12.2.1/python --cpp-dir tests/benchmarks/baselines/cpu/amd-epyc-9534_96/gcc_12.2.1/cpp
```

## Compare results

Compare current cpp benchmarks against baseline
```bash
./scripts/benchmarks/compare_benchmarks.sh tests/benchmarks/baselines/cpu/amd-epyc-9534_96/gcc_12.2.1/cpp build/Benchmarking/cpp
```

