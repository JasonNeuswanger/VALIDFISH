//
// Created by Jason Neuswanger on 2/21/18.
//

#include "third-party/pybind11/include/pybind11/pybind11.h"
#include "third-party/pybind11/include/pybind11/stl.h"
#include "Forager.h"
#include "Optimizer.h"

namespace py = pybind11;

PYBIND11_MODULE(pyvalidfish, m) {

    m.doc() = "Visual Attention-Limited Individual Drift Forager In Simulated Habitat (VALIDFISH) model in C++."; // optional module docstring

    /* --------------------- Prey category class and components ----------------------- */

    py::class_<PreyCategory> preycat(m, "PreyCategory");
    preycat.def(py::init<int, std::string, double, double, double, double, double, double>())
            .def("get_diet_proportion", &PreyCategory::get_diet_proportion);

    /* --------------------- Forager class and components ----------------------- */

    py::class_<Forager> forager(m, "Forager");
    forager.def(py::init<double, double, double, double, double, double, double, double, double, double, double, double,
                    double, double, unsigned, double, double, double, double, std::string*>())
            .def("NREI", &Forager::NREI)
            .def("GREI", &Forager::GREI)
            .def("modify_strategy", &Forager::modify_strategy)
            .def("modify_parameter", &Forager::modify_parameter)
            .def("modify_parameters", &Forager::modify_parameters)
            .def("modify_bound", &Forager::modify_bound)
            .def("fix_bound", &Forager::fix_bound)
            .def("add_prey_category", &Forager::add_prey_category)
            .def("process_prey_category_changes", &Forager::process_prey_category_changes)
            .def("depletable_detection_probability", &Forager::depletable_detection_probability)
            .def("relative_pursuits_by_position", &Forager::relative_pursuits_by_position)
            .def("relative_pursuits_by_position_single_prey_type", &Forager::relative_pursuits_by_position_single_prey_type)
            .def("proportion_of_detections_within", &Forager::proportion_of_detections_within)
            .def("get_prey_category", &Forager::get_prey_category)
            .def("get_fork_length_cm", &Forager::get_fork_length_cm)
            .def("get_radius", &Forager::get_radius)
            .def("get_theta", &Forager::get_theta)
            .def("get_focal_velocity", &Forager::get_focal_velocity)
            .def("get_foraging_attempt_rate", &Forager::get_foraging_attempt_rate)
            .def("get_proportion_of_attempts_ingested", &Forager::get_proportion_of_attempts_ingested)
            .def("get_diet_proportion_for_prey_category", &Forager::get_diet_proportion_for_prey_category)
            .def("get_angular_resolution", &Forager::get_angular_resolution)
            .def("print_strategy", &Forager::print_strategy)
            .def("print_status", &Forager::print_status)
            .def("print_analytics", &Forager::print_analytics)
            .def("analyze_results", &Forager::analyze_results);

    py::enum_<Forager::Strategy>(forager, "Strategy")
            .value("radius", Forager::Strategy::s_radius)
            .value("theta", Forager::Strategy::s_theta)
            .value("mean_column_velocity", Forager::Strategy::s_mean_column_velocity)
            .value("saccade_time", Forager::Strategy::s_saccade_time)
            .value("discrimination_threshold", Forager::Strategy::s_discrimination_threshold)
            .export_values();

    py::enum_<Forager::Parameter>(forager, "Parameter")
            .value("delta_0", Forager::Parameter::p_delta_0)
            .value("alpha_0", Forager::Parameter::p_alpha_0)
            .value("beta", Forager::Parameter::p_beta)
            .value("Z_0", Forager::Parameter::p_Z_0)
            .value("c_1", Forager::Parameter::p_c_1)
            .export_values();

    /* --------------------- Optimizer class and components ----------------------- */

    py::class_<Optimizer> optimizer(m, "Optimizer");
    optimizer.def(py::init<Forager *, size_t, size_t, bool>())
            .def("optimize_forager", &Optimizer::optimize_forager)
            .def("set_algorithm_options", &Optimizer::set_algorithm_options);

}
