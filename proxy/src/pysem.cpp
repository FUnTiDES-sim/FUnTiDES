#include <nanobind/nanobind.h>
#include "SEMproxy.hpp"

int TOTO(int a, int b) { return a + b; }

NB_MODULE(pysem, m) {
    m.def("TOTO", &TOTO);
}
