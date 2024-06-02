#include <pybind11/pybind11.h>

namespace py = pybind11;

extern void init_Nuclide(py::module &);

extern void init_Chain(py::module &);

PYBIND11_MODULE(hericendre, m) {
    init_Nuclide(m);
    init_Chain(m);
}
