import pysem

pysem.initialize_kokkos()

sem = pysem.SEMproxy(30,30,30,1000)
sem.initFiniteElem()
sem.run()

