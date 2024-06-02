#include <chain.h>
#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

void init_Chain(py::module &m) {
    m.doc() = "Chain class";

    py::class_<Chain>(m, "Chain").def(py::init<std::string>());
}
