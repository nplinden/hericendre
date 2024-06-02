#include <nuclide.h>
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

void init_Nuclide(py::module &m) {
    m.doc() = "Nuclide class";

    py::class_<Nuclide>(m, "Nuclide")
            .def(py::init<std::string, double>())
            .def_property("name", &Nuclide::getName, &Nuclide::setName)
            .def_property("dconst", &Nuclide::getDconst, &Nuclide::setDconst);
}
