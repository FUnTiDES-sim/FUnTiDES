#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import kokkos
import time
from datetime import datetime

import pysolver as Solver


def print_global(arr):
    for idx, e in enumerate(arr):
        if e != 0:
            print(f"arr[{idx}] = {e}")


def sourceTerm(time_n, f0):
    o_tpeak = 1.0 / f0
    pulse = 0.0
    if time_n <= -0.9 * o_tpeak or time_n >= 2.9 * o_tpeak:
        return pulse

    pi = 3.14157
    lam = (f0 * pi) * (f0 * pi)
    pulse = (
        2.0
        * lam
        * (2.0 * lam * (time_n - o_tpeak) * (time_n - o_tpeak) - 1.0)
        * np.exp(-lam * (time_n - o_tpeak) * (time_n - o_tpeak))
    )
    return pulse


def initModel(model, nElements):
    for i in range(nElements):
        model[i] = 1500.0
    return 0


def initPressure(pressure, nDof):
    for i in range(nDof):
        pressure[i, 0] = 0.0
        pressure[i, 1] = 0.0
    return 0


def getSnapshot(timeStep, i1, nx, ny, nz, hx, hy, hz, pnGlobal, normalize):
    offset = nx * nz * (int(ny / 2) - 1)
    grid = np.zeros((nx, nz))
    for I in range(offset, offset + nx * nz):
        i = (I - offset) % nx
        j = int((I - offset - i) / nx)
        grid[i, j] = pnGlobal[I, i1]

    if normalize:
        maxvalue = np.abs(grid).max()
        if maxvalue != 0:
            grid = grid / maxvalue

    return grid


def get_element_number_from_point(ex, ey, ez, hx, hy, hz, x, y, z):
    eX = eY = eZ = 0
    for i in range(ex):
        if x >= i * hx and x < (i + 1) * hx:
            eX = i
    for i in range(ey):
        if y >= i * hy and y < (i + 1) * hy:
            eY = i
    for i in range(ez):
        if z >= i * hz and z < (i + 1) * hz:
            eZ = i

    return eX + eZ * ex + eY * ex * ez


def main():
    kokkos.initialize()
    # Add timing variables
    start_time = time.time()
    simulation_start = datetime.now()
    iteration_times = []

    print(f"Simulation started at: {simulation_start}")

    order = 2
    domain_size = 1500.0
    ex = 200
    ey = 200
    ez = 200
    hx = domain_size / ex
    hy = domain_size / ey
    hz = domain_size / ez

    nx = ex * order + 1
    ny = ey * order + 1
    nz = ez * order + 1

    nDof = nx * ny * nz

    nElements = ex * ey * ez
    nPointsPerElement = (order + 1) * (order + 1) * (order + 1)
    freqPrint = 5
    printOffset = int(nDof / freqPrint)
    print("number of elements=", nElements)

    memspace = kokkos.CudaUVMSpace
    layout = kokkos.LayoutLeft

    kk_model = kokkos.array(
        [nElements], dtype=kokkos.float32, space=memspace, layout=layout
    )
    model = np.array(kk_model, copy=False)
    initModel(model, nElements)

    mesh = Solver.SEMmesh(ex, ey, ez, domain_size, domain_size, domain_size, order)

    solver = Solver.SEMsolver(mesh)

    nodesList = np.zeros((nElements, nPointsPerElement), dtype=int)

    # allocate pressure
    print("Allocating Pressure...")
    kk_pnGlobal = kokkos.array(
        [nDof, 2], dtype=kokkos.float32, space=memspace, layout=layout
    )
    pnGlobal = np.array(kk_pnGlobal, copy=False)
    initPressure(pnGlobal, nDof)

    # source term
    f0 = 5
    # time step and sampling
    timeStep = 0.001
    nTimeSteps = 3000
    numberOfRHS = 1
    xs = ex * hx / 2
    ys = ey * hy / 2
    zs = ez * hz / 2

    kk_RHSElement = kokkos.array(
        [numberOfRHS], dtype=kokkos.int32, space=memspace, layout=layout
    )
    RHSElement = np.array(kk_RHSElement, copy=False)
    RHSElement[0] = get_element_number_from_point(ex, ey, ez, hx, hy, hz, xs, ys, zs)
    print("RHS element number ", RHSElement[0])

    # compute source term
    kk_RHSTerm = kokkos.array(
        [numberOfRHS, nTimeSteps], dtype=kokkos.float32, space=memspace, layout=layout
    )
    RHSTerm = np.array(kk_RHSTerm, copy=False)
    for i in range(nTimeSteps):
        RHSTerm[0, i] = sourceTerm(i * timeStep, f0)

    # setup graphic display
    grid = np.zeros((nx, nz))
    fig, ax = plt.subplots()
    cmpvalue = 20
    im = ax.imshow(
        grid, cmap="viridis", interpolation="nearest", vmin=-cmpvalue, vmax=cmpvalue
    )

    plt.colorbar(im, ax=ax, label="Intensity")
    plt.title("2D Slice of a Float32 Array")
    plt.xlabel("X-axis")
    plt.ylabel("Z-axis")
    # propagate one step forward in time
    i1 = 0
    i2 = 1
    maxval = 1

    for timeSample in range(nTimeSteps):
        iter_start = time.time()
        solver.computeOneStep(
            timeSample,
            i1,
            i2,
            kk_RHSTerm,
            kk_pnGlobal,
            kk_RHSElement,
        )

        iter_time = time.time() - iter_start
        iteration_times.append(iter_time)

        if timeSample % 1000 == 0:
            print("sum pnGlobal[:, i1]", np.sum(pnGlobal[:, i1]))
            print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
            elementSource = nodesList[RHSElement[0], 0]
            elapsed_time = time.time() - start_time
            print()

        if timeSample % 100 == 0:
            print(f"Time {timeSample} / {nTimeSteps}")

        if timeSample % 10 == 0:
            # plotting
            grid = getSnapshot(timeStep, i1, nx, ny, nz, hx, hy, hz, pnGlobal, False)
            im.set_array(grid)  # Update plot with new values
            plt.draw()  # Redraw the figure with updated data
            plt.ioff
            plt.savefig(f"snap0{timeSample:0{5}d}.png")

        # Swap pn and pn+1
        tmp = i1
        i1 = i2
        i2 = tmp

    # Print final timing statistics
    end_time = time.time()
    simulation_end = datetime.now()
    total_time = end_time - start_time

    print("\nSimulation Statistics:")
    print(f"Start time: {simulation_start}")
    print(f"End time: {simulation_end}")
    print(f"Total runtime: {total_time:.2f} seconds")
    print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
    print(f"Min iteration time: {np.min(iteration_times):.4f} seconds")
    print(f"Max iteration time: {np.max(iteration_times):.4f} seconds")

    # release kokkos arrays and vectors
    del kk_pnGlobal
    del kk_model
    del kk_RHSTerm
    del kk_RHSElement
    del solver
    del mesh

    kokkos.finalize()
    print("end of  computation")
    print("we close kokkos")


if __name__ == "__main__":
    main()
