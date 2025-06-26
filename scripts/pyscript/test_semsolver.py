#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import libpykokkos as kokkos
import time
from datetime import datetime

import pysolver as Solver
import pysem as Sem

ArrayReal = kokkos.KokkosView_float32_HostSpace_LayoutRight_2
VectorInt = kokkos.KokkosView_int32_HostSpace_LayoutRight_1
VectorReal = kokkos.KokkosView_float32_HostSpace_LayoutRight_1

# ArrayReal = kokkos.KokkosView_float32_CudaUVMSpace_LayoutLeft_2
# VectorInt = kokkos.KokkosView_int32_CudaUVMSpace_LayoutLeft_1
# VectorReal = kokkos.KokkosView_float32_CudaUVMSpace_LayoutLeft_1

def print_global( arr ):
  for idx, e in enumerate(arr):
    if e != 0:
      print(f"arr[{idx}] = {e}")

def sourceTerm( time_n, f0 ):
  o_tpeak = 1.0/f0;
  pulse = 0.0;
  if time_n <= -0.9*o_tpeak or time_n >= 2.9*o_tpeak:
    return pulse

  pi = 3.14157
  lam = (f0*pi)*(f0*pi)
  pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*np.exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak))
  return pulse

def initModel(model,nElements):
    for i in range(nElements):
        model[i]=1500.
    return 0

def initPressure(pressure,nDof):
    for i in range(nDof):
        pressure[i,0]=0.
        pressure[i,1]=0.
    return 0

def getSnapshot(timeStep,i1,nx,ny,nz,hx,hy,hz,pnGlobal,maxval):
  offset=nx*nz*(int(ny/2)-1)
  print(nx,ny,nz,offset)
  grid=np.zeros((nx,nz))
  for I in range(offset,offset+nx*nz):
      i=(I-offset)%nx
      j=int((I-offset-i)/nx);
      if maxval != 0:
        grid[i,j]=pnGlobal[I,i1]/maxval
      else :
        grid[i,j]=pnGlobal[I,i1]
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

  # kokkos.initialize()

  order=2
  domain_size=1500.
  ex=100
  ey=100
  ez=100
  hx=domain_size/ex
  hy=domain_size/ey
  hz=domain_size/ez

  nx=ex*order+1
  ny=ey*order+1
  nz=ez*order+1

  nDof=nx*ny*nz

  nElements=ex*ey*ez
  nPointsPerElement=(order+1)*(order+1)*(order+1)
  freqPrint=5
  printOffset=int(nDof/freqPrint)
  print("number of elements=",nElements)

  model = np.zeros([nElements],dtype=np.float32)
  initModel(model,nElements)
  kk_model = VectorReal(model, (nElements,))

  mesh = Sem.SEMmesh(ex, ey, ez, domain_size, domain_size, domain_size, order, 20, False)
  myInfo = Sem.SEMinfo()
  myInfo.numberOfNodes = nx * ny * nz
  myInfo.numberOfElements = nElements
  myInfo.numberOfPointsPerElement = nPointsPerElement
  myInfo.numberOfInteriorNodes = mesh.getNumberOfInteriorNodes()

  solver = Solver.SEMsolver()
  solver.computeFEInit(myInfo, mesh)
  solver.set_model(kk_model)

  nodesList=np.zeros((nElements,nPointsPerElement),dtype=int)

  # allocate pressure
  print("Allocating Pressure...")
  pnGlobal = np.zeros((nDof, 2), dtype=np.float32)
  initPressure(pnGlobal, nDof)
  kk_pnGlobal = ArrayReal(pnGlobal, (nDof, 2))

  # source term
  f0=5
  # time step and sampling
  timeStep=0.001
  nTimeSteps=300
  numberOfRHS=1
  xs=ex*hx/2
  ys=ey*hy/2
  zs=ez*hz/2

  RHSElement=np.array([numberOfRHS],dtype=int)
  RHSElement[0] = get_element_number_from_point(ex,ey,ez,hx,hy,hz,xs,ys,zs)
  kk_RHSElement = VectorInt(RHSElement, (numberOfRHS,))
  print("RHS element number ", RHSElement[0])

  # compute source term
  RHSTerm=np.zeros((numberOfRHS,nTimeSteps),dtype=np.float32)
  for i in range(nTimeSteps):
      RHSTerm[0,i]=sourceTerm(i*timeStep,f0)
  kk_RHSTerm = ArrayReal(RHSTerm, (numberOfRHS,nTimeSteps))

  # setup graphic display
  grid=np.zeros((nx,nz))
  plt.ion()  # Enable interactive mode
  fig, ax = plt.subplots()
  im = ax.imshow(grid, cmap='viridis', interpolation='nearest')
  plt.colorbar(im, ax=ax, label="Intensity")
  plt.title("2D Slice of a Float32 Array")
  plt.xlabel("X-axis")
  plt.ylabel("Z-axis")
  # propagate one step forward in time
  i1=0
  i2=1
  maxval=1

  for timeSample in range(nTimeSteps):
     iter_start = time.time()
     if timeSample%100==0:
       print()
       print("sum pnGlobal[:, i1]", np.sum(pnGlobal[:, i1]))

     solver.computeOneStep(timeSample, order, nPointsPerElement, i1, i2, myInfo, kk_RHSTerm, kk_pnGlobal, kk_RHSElement)
     # Solver.compute_debug_args(solver, 0, 2, 27, 0, 1, myInfo,kk_RHSTerm, kk_pnGlobal, kk_RHSElement)

     iter_time = time.time() - iter_start
     iteration_times.append(iter_time)

     if timeSample%10==0:
       print("sum pnGlobal[:, i1]", np.sum(pnGlobal[:, i1]))
       maxval=np.max(np.abs(pnGlobal))/5
       print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
       elementSource=nodesList[RHSElement[0],0]
       elapsed_time = time.time() - start_time
       print("TIME")
       print(f"Time iteration {timeSample}/{nTimeSteps}")
       print(f"Elapsed time: {elapsed_time:.2f} seconds")
       print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")
       print(f"Pressure={pnGlobal[elementSource,0]}")

       grid=getSnapshot(timeStep,i1,nx,ny,nz,hx,hy,hz,pnGlobal,maxval)
       im.set_array(grid)  # Update plot with new values
       plt.draw()  # Redraw the figure with updated data
       plt.pause(1)  # Pause for 1 second before updating again
       plt.ioff
       plt.show()

     tmp=i1
     i1=i2
     i2=tmp

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
  # del nodesList, nodesCoordsX, nodesCoordsY,nodesCoordsZ,model,pnGlobal
  # del massMatrixGlobal,yGlobal,RHSElement,RHSTerm
  del solver
  del mesh
  del myInfo

  kokkos.finalize()
  print( "end of  computation")
  print("we close kokkos")

if __name__ == "__main__":
  main()
