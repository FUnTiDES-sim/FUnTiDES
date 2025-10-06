#ifndef FDTDMACROS_HPP
#define FDTDMACROS_HPP

#define NX myGrids.nx
#define NY myGrids.ny
#define NZ myGrids.nz

#define LX myGrids.lx
#define LY myGrids.ly
#define LZ myGrids.lz

#define POW2(x) ((x) * (x))
#define IDX3(i, j, k) (NZ * NY * (i) + NZ * (j) + (k))
#define IDX3_l(i, j, k)                                                      \
  ((NZ + 2 * LZ) * (NY + 2 * LY) * ((i) + LX) + (NZ + 2 * LZ) * ((j) + LY) + \
   ((k) + LZ))
#define IDX3_eta1(i, j, k) \
  ((NZ + 2) * (NY + 2) * ((i) + 1) + (NZ + 2) * ((j) + 1) + ((k) + 1))

#if defined(USE_KOKKOS)
#define LOOP3DHEAD(x3, y3, z3, x4, y4, z4) Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<3>>({z3,x3,y3},{z4,x4,y4}),KOKKOS_LAMBDA(int k,int i,int j) {
#else
#define LOOP3DHEAD(x3, y3, z3, x4, y4, z4) \
  for (int i = x3; i < x4; ++i)            \
  {                                        \
    for (int j = y3; j < y4; ++j)          \
    {                                      \
      for (int k = z3; k < z4; ++k)        \
      {
#endif

#if defined(USE_KOKKOS)
#define LOOP3DEND \
  });
#else
#define LOOP3DEND \
  }               \
  }               \
  }
#endif

#define CREATEVIEWINNER              \
  vectorReal coefx = myModels.coefx; \
  vectorReal coefy = myModels.coefy; \
  vectorReal coefz = myModels.coefz; \
  vectorReal vp = myModels.vp;       \
  double coef0 = myModels.coef0;
#define CREATEVIEWPML            \
  CREATEVIEWINNER                \
  vectorReal phi = myModels.phi; \
  vectorReal eta = myModels.eta;
#define CREATEVIEWRHS                    \
  vectorReal RHSTerm = myModels.RHSTerm; \
  vectorReal pn = myModels.pn;           \
  vectorReal vp = myModels.vp;
#define CREATEVIEWSPONGE vectorReal spongeArray = myModels.spongeArray;
#define PN_Global pnGlobal

#if defined(USE_KOKKOS)
#define FDFENCE Kokkos::fence();
#else
#define FDFENCE
#endif

#endif  // FDTDMACROS_HPP
