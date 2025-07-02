#!/usr/bin/env python3
import sys
import faulthandler
faulthandler.enable()
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

#ArrayReal = kokkos.KokkosView_float32_HostSpace_LayoutLeft_2
#VectorInt = kokkos.KokkosView_int32_HostSpace_LayoutLeft_1
#VectorReal = kokkos.KokkosView_float32_HostSpace_LayoutLeft_1

def print_global( arr ):
  for idx, e in enumerate(arr):
    if e != 0:
      print(f"arr[{idx}] = {e}")

def sourceTerm(time_n, f0):
    o_tpeak = 1.2/f0 # golden rule (it just works this way)
    # This is close to 2.0*np.sqrt(np.pi) / (3*f0), where a cut-off frequency of 3*f0 is assumed
    pulse = 0.0

    # Standard truncation for Ricker wavelet
    if time_n <= -0.9*o_tpeak or time_n >= 2.9*o_tpeak:
        return pulse

    pi = 3.14157
    lam = (f0*pi)*(f0*pi)
    pulse = 2.0*lam*(2.0*lam*(time_n-o_tpeak)*(time_n-o_tpeak)-1.0)*np.exp( -lam*(time_n-o_tpeak)*(time_n-o_tpeak))
    return pulse

def initModel(model,nElements):
    for i in range(nElements):
        model[i]=2000.
    return 0

def initPressure(pressure,nDof):
    for i in range(nDof):
        pressure[i,0]=0.
        pressure[i,1]=0.
    return 0

def getSnapshot(i1,xs,I,timeSample,ny,nz,pnGlobal,maxval):

  pressure_vals = pnGlobal[I, i1].copy()

  # Normalize pressure values for better visualization
  if maxval  > 0:
    pressure_vals = pressure_vals / maxval

  #pressure_vals = np.reshape(pressure_vals, (124*self.order+1, 123*self.order+1))
  pressure_vals = np.reshape(pressure_vals, (ny,nz))

  maxval = np.max(np.abs(pressure_vals))/50
  pressure_vals = np.clip(pressure_vals, -maxval, maxval)

  plt.figure(figsize=(8, 6))
  plt.imshow(pressure_vals.T , origin='lower', cmap='viridis', interpolation='nearest')
  plt.colorbar(label='Pressure')
  plt.xlabel('Y')
  plt.ylabel('Z')
  plt.title(f'Pressure field at X={xs:.3f}, timeSample={timeSample}')
  plt.tight_layout()
  file_name = f"semProxyMESH_snapshot_{timeSample}.png"
  plt.savefig(file_name, dpi=150)
  plt.close()
  print(f"Snapshot PNG file written as {file_name}", flush=True)
  

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

  print("Propagating...", flush=True)

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

  print("number of elements=",nElements)

  model = np.zeros([nElements],dtype=np.float32)
  initModel(model,nElements)
  kk_model = VectorReal(model, (nElements,))

  mesh = Sem.SEMmesh(ex, ey, ez, ex*hx, ey*hy, ez*hz, order, 0, False)
  myInfo = Sem.SEMinfo()
  myInfo.numberOfNodes = nx * ny * nz
  myInfo.numberOfElements = nElements
  myInfo.numberOfPointsPerElement = nPointsPerElement
  myInfo.numberOfInteriorNodes = mesh.getNumberOfInteriorNodes()
  myInfo.numberOfSpongeNodes = mesh.getSpongeSize()
  
  print(dir(myInfo))

  solver = Solver.SEMsolver()
  solver.computeFEInit(myInfo, mesh)
  solver.set_model(kk_model)

  # allocate pressure
  print("Allocating Pressure...")
  pnGlobal = np.zeros((nDof, 2), dtype=np.float32)
  initPressure(pnGlobal, nDof)
  kk_pnGlobal = ArrayReal(pnGlobal, (nDof, 2))

  # source term
  f0=5 # Peak frequency of a Ricker wavelet
  # time step and sampling
  timeStep=0.001
  nTimeSteps=300
  numberOfRHS=1
  xs=735
  ys=735
  zs=705

  maxAmplitude = 1000

  RHSElement=np.array([numberOfRHS],dtype=np.int32)
  RHSElement[0] = 474949 #get_element_number_from_point(ex,ey,ez,hx,hy,hz,xs,ys,zs)
  kk_RHSElement = VectorInt(RHSElement, (numberOfRHS,))
  print("RHS element number ", RHSElement[0])

  # compute source term
  RHSTerm = np.zeros((numberOfRHS, nTimeSteps), dtype=np.float32)
  for i in range(nTimeSteps):
      RHSTerm[0, i] = sourceTerm(i*timeStep, f0)
  RHSTerm = RHSTerm/np.max(np.abs(RHSTerm))
  RHSTerm *= maxAmplitude
  kk_RHSTerm = ArrayReal(RHSTerm, (numberOfRHS, nTimeSteps))

  nodesList = np.empty((nElements, nPointsPerElement), dtype=np.int32)
  xi = [-1,0,1]
  nodesCoordsX = np.empty((nx * ny * nz,), dtype=np.float32)
  
  for j in range(ey):
      for k in range(ez):
          for i in range(ex):
              n0 = i + k * ex + j * ex * ez
              offset = i * order + k * order * nx + j * order * nx * nz

              x0 = i * hx
              x1 = (i + 1) * hx
              b = (x1 + x0) / 2.
              a = b - x0

              for m in range(order + 1):
                  for n in range(order + 1):
                      for l in range(order + 1):
                          dofLocal = l + n * (order + 1) + m * (order + 1) * (order + 1)
                          dofGlobal = offset + l + n * nx + m * nx * nz
                          nodesList[n0, dofLocal] = dofGlobal
                          nodesCoordsX[dofGlobal] = a * xi[l] + b


  print(f"dofGlobal : {dofGlobal}", flush=True)

  rhs_node = nodesList[RHSElement[0], 0]
  print(f"RHS node number {rhs_node}", flush=True)

  tol = 1e-6
  mask = np.abs(nodesCoordsX - xs) < tol
  plane_4_snapshot = np.where(mask)[0]

  print(f"size of nodesCoordsX: {nodesCoordsX.shape}", flush=True)
  print(f"nodesCoordsX: {nodesCoordsX}", flush=True)
  print(f"shape of plane_4_snapshot: {plane_4_snapshot.shape}", flush=True)
  print(f"Plane 4 snapshot indices: {plane_4_snapshot}", flush=True)

  #sys.exit(0)

  # Add timing variables
  start_time = time.time()
  simulation_start = datetime.now()
  iteration_times = []

  print(f"Simulation started at: {simulation_start}", flush=True)

  # propagate one step forward in time
  i1=0
  i2=1
  maxval=1

  for timeSample in range(nTimeSteps):
     print(f"   - time iteration {timeSample}/{nTimeSteps}", flush=True)
     print(f"i1 = {i1}, i2 = {i2}", flush=True)
     iter_start = time.time()

     solver.computeOneStep(timeSample, order, nPointsPerElement, i1, i2, myInfo, kk_RHSTerm, kk_pnGlobal, kk_RHSElement)
     
     iter_time = time.time() - iter_start
     iteration_times.append(iter_time)

     if timeSample % 10 == 0:
        print(f"     sum pnGlobal[:, i1] = {np.sum(pnGlobal[:, i1])}")
        print(f"     Percentage of zeros in pnGlobal = {np.count_nonzero(pnGlobal[:,i1] == 0) / pnGlobal[:,i1].size * 100}", flush=True)
        print(f"     RHSElement[0] = {RHSElement[0]}")
        maxval = np.max(np.abs(pnGlobal[:,i1]))              
        getSnapshot(i1,xs,plane_4_snapshot,timeSample,ny,nz,pnGlobal,maxval)
        print(f"     Pressure max val after iteration {timeSample}: {maxval}")
        print(f"     Average iteration time: {np.mean(iteration_times):.4f} seconds")
        elapsed_time = time.time() - start_time
        print(f"     Time iteration {timeSample}/{nTimeSteps}")
        print(f"     Elapsed time: {elapsed_time:.2f} seconds")
        print(f"     Average iteration time: {np.mean(iteration_times):.4f} seconds", flush=True)

     i1, i2 = i2, i1

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
  print(f"Max iteration time: {np.max(iteration_times):.4f} seconds", flush=True)

  # release kokkos arrays and vectors
  # del nodesList, nodesCoordsX, nodesCoordsY,nodesCoordsZ,model,pnGlobal
  # del massMatrixGlobal,yGlobal,RHSElement,RHSTerm
  del solver
  del mesh
  del myInfo

  kokkos.finalize()
  print("THE END")
  print("THAT IS ALL FOLKS!")

if __name__ == "__main__":
  main()
