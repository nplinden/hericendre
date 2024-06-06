#include <model.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>

namespace py = pybind11;

void init_Model(py::module &m) {
    m.doc() = "Model class";

    py::class_<Model>(m, "Model")
            .def(py::init<std::string>())
            .def(py::init<>())
            .def("run", &Model::run)
            .def_readwrite("resultpath", &Model::resultpath_)
            .def_readwrite("solvertype", &Model::solvertype_)
            .def_readwrite("times", &Model::times_)
            .def_readwrite("initcc", &Model::initcc_)
            .def_property("chainpath", &Model::chainpath, &Model::set_chainpath);
}
