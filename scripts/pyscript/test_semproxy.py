import pysem
import libpykokkos as kokkos

# pysem.initialize_kokkos()
kokkos.initialize()

sem = pysem.SEMproxy(100,100,100,2000)
sem.initFiniteElem()
sem.run()
del sem
# pysem.finalize_kokkos()
kokkos.finalize()
