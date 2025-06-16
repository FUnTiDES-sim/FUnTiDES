#!/usr/bin/env python3
import faulthandler
faulthandler.enable()
import numpy as np
import matplotlib.pyplot as plt
# import libpykokkos as kokkos
# import kokkos
import time
from datetime import datetime

import pysolver as Solver
import pysem as Sem

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
      grid[i,j]=pnGlobal[I,i1]/maxval
  return grid

def main():

  # Add timing variables
  start_time = time.time()
  simulation_start = datetime.now()
  iteration_times = []

  print(f"Simulation started at: {simulation_start}")

  # kokkos.initialize()

  order=2
  ex=100
  ey=100
  ez=100
  hx=15.
  hy=15.
  hz=15.

  nx=ex*order+1
  ny=ey*order+1
  nz=ez*order+1

  nDof=nx*ny*nz

  nElements=ex*ey*ez
  nPointsPerElement=(order+1)*(order+1)*(order+1)
  freqPrint=5
  printOffset=int(nDof/freqPrint)
  print("number of elements=",nElements)

  mesh = Sem.SEMmesh(ex, ey, ez, hx, hy, hz, order, 0, False)
  myInfo = Sem.SEMinfo()
  myInfo.numberOfNodes = nx * ny * nz
  myInfo.numberOfElements = nElements
  myInfo.numberOfPointsPerElement = nPointsPerElement

  solver = Solver.SEMsolver()
  solver.initFEarrays(myInfo, mesh)

  nodesList=np.zeros((nElements,nPointsPerElement),dtype=int)
  # SEM.getGlobalNodesList(order,ex,ey,ez,nx,ny,nz,nodesList)

  # nodesCoordsX=kokkos.array([nElements,nPointsPerElement],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # nodesCoordsY=kokkos.array([nElements,nPointsPerElement],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # nodesCoordsZ=kokkos.array([nElements,nPointsPerElement],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # SEM.getNodesCoordinates(order,ex,ey,ez,hx,hy,hz,nodesCoordsX,nodesCoordsY,nodesCoordsZ)

  # allocate model
  model = np.zeros(nElements, dtype=np.float32)
  initModel(model, nElements)

  # allocate pressure
  pnGlobal = np.zeros((nDof, 2), dtype=np.float32)
  initPressure(pnGlobal, nDof)

  # allocate mass matrix and stiffness vector
  massMatrixGlobal = np.zeros(nDof, dtype=np.float32)
  yGlobal = np.zeros(nDof, dtype=np.float32)

  # # allocate model
  # model=kokkos.array([nElements],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # initModel(model,nElements)

  # # allocate pressure
  # pnGlobal=kokkos.array([nDof,2],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # initPressure(pnGlobal,nDof)

  # # allocate mass matrix and stiffness vector
  # massMatrixGlobal=kokkos.array([nDof],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # yGlobal=kokkos.array([nDof],dtype=kokkos.float,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)

  # source term
  f0=15
  # time step and sampling
  timeStep=0.001
  timeStep2=timeStep*timeStep
  nTimeSteps=1000
  numberOfRHS=1
  xs=ex*hx/2
  ys=ey*hy/2
  zs=ez*hz/2

  RHSElement=np.array([numberOfRHS],dtype=int)
  # RHSElement=kokkos.array([numberOfRHS],dtype=kokkos.int,layout=kokkos.LayoutLeft,space=kokkos.CudaUVMSpace)
  # RHSElement[0]=SEM.getElementNumberFromPoint(ex,ey,ez,hx,hy,hz,xs,ys,zs)
  print("RHS elemet number ",RHSElement[0])

  # compute source term
  RHSTerm=np.zeros((numberOfRHS,nTimeSteps),dtype=float)
  for i in range(nTimeSteps):
      RHSTerm[0,i]=sourceTerm(i*timeStep,f0)

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

     # SEM.computeOneStep(timeSample, order, nPointsPerElement, i1, i2, myInfo, RHSTerm, pnGlobal, RHSElement)

     iter_time = time.time() - iter_start
     iteration_times.append(iter_time)

     if timeSample==100:
        maxval=np.max(np.abs(pnGlobal))/5
        print(f"Pressure max val after iteration 100: {maxval}")
        print(f"Average iteration time: {np.mean(iteration_times):.4f} seconds")

     if timeSample%100==0:
        elementSource=nodesList[RHSElement[0],0]
        elapsed_time = time.time() - start_time
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
  del nodesList, nodesCoordsX, nodesCoordsY,nodesCoordsZ,model,pnGlobal
  del massMatrixGlobal,yGlobal,RHSElement,RHSTerm

  # kokkos.finalize()
  print( "end of  computation")
  print("we close kokkos")

if __name__ == "__main__":
  main()
