#include <model.h>
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

void init_Model(py::module &m) {
    m.doc() = "Model class";

    py::class_<Model>(m, "Model").def(py::init<std::string>());
}
