#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "protein2.h"
#include "sp_type.h"


namespace py = pybind11;


void get_u(Salign& salign, py::array_t<double, py::array::c_style>& array){
    auto ptr = static_cast<double*>(array.request().ptr);
    //u[4][3]

    double (*u)[4][3] = salign.get_u();

    for (int i=0; i<4; i++){
        for (int j=0; j<3; j++){
            ptr[i*3+j] = (*u)[i][j];
        }
    }
}

PYBIND11_MODULE(SPlib, m) {
    m.doc() = "SPlib";
    py::class_<Protein2>(m, "Protein2")
        .def(py::init<int, std::string &, vector<vector<double> > >())
        .def(py::init<std::string &>())
        .def("calSS", &Protein2::calSS)
        .def("calSS2", &Protein2::calSS2)
        .def("add_ideal", &Protein2::add_ideal)
        .def("add_dummy_ideal", &Protein2::add_dummy_ideal)
        .def("getsegstart", &Protein2::getsegstart)
        .def("getsegmid", &Protein2::getsegmid)
        .def("getsegend", &Protein2::getsegend)
        .def("getssec", &Protein2::getssec);
    py::class_<Salign>(m, "Salign")
        .def(py::init<>())
        .def("Run_all", &Salign::Run_all)
        .def("RunSeg", &Salign::RunSeg)
        .def("init", static_cast<void (Salign::*)(Protein2 *, Protein2 *)>(&Salign::init))
        .def("calscores_all", &Salign::calscores_all)
        .def("getscore", &Salign::getscore)
        .def("calRscore", &Salign::calRscore)
        .def("getscoresall", &Salign::getscoresall)
        .def("score_final", &Salign::score_final)
        .def("restore_max", &Salign::restore_max)
        .def("calSPscore", &Salign::calSPscore)
        .def("calSPscoreQ", &Salign::calSPscoreQ)
        .def("run_pairwise", static_cast<double (Salign::*)(Protein2 *, Protein2 *)>(&Salign::run_pairwise))
        .def("getRmax", &Salign::getRmax);
    m.def("get_u", &get_u);
}


