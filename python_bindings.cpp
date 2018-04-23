//
// Created by Jason Neuswanger on 2/21/18.
//

#include "lib/pybind11/include/pybind11/pybind11.h"
#include "lib/pybind11/include/pybind11/stl.h"
#include "Forager.h"
#include "Optimizer.h"

namespace py = pybind11;

PYBIND11_MODULE(pyvalidfish, m) {

    m.doc() = "Visual Attention-Limited Individual Drift Forager In Simulated Habitat (VALIDFISH) model in C++."; // optional module docstring

    /* --------------------- Prey type class and components ----------------------- */

    py::class_<PreyType> preytype(m, "PreyType");

    preytype.def(py::init<int, std::string, double, double, double, double, bool, double>())
            .def("get_diet_proportion", &PreyType::get_diet_proportion)
            .def("get_false_positive_probability", &PreyType::get_false_positive_probability)
            .def("get_true_hit_probability", &PreyType::get_true_hit_probability)
            .def("get_perceptual_sigma", &PreyType::get_perceptual_sigma)
            .def("get_max_visible_distance", &PreyType::get_max_visible_distance)
            .def("get_max_attended_distance", &PreyType::get_max_attended_distance)
            .def("get_name", &PreyType::get_name)
            .def("get_length", &PreyType::get_length);

    /* --------------------- Forager class and components ----------------------- */

    py::class_<Forager> forager(m, "Forager");

    forager.def(py::init<double, double, double, double, double, double, double, double, double, double, double, double,
            double, double, double, double, unsigned, double, double , double, double, double, std::string *>())
            .def("NREI", &Forager::NREI)
            .def("GREI", &Forager::GREI)
            .def("tau", &Forager::tau)
            .def("tau_components", &Forager::tau_components)
            .def("detection_probability", &Forager::detection_probability)
            .def("passage_time", &Forager::passage_time)
            .def("set_strategy", &Forager::set_strategy)
            .def("set_strategies", &Forager::set_strategies)
            .def("set_parameter", &Forager::set_parameter)
            .def("set_parameters", &Forager::set_parameters)
            .def("set_single_strategy_bounds", &Forager::set_single_strategy_bounds)
            .def("fix_single_strategy_bound", &Forager::fix_single_strategy_bound)
            .def("add_prey_type", &Forager::add_prey_type)
            .def("process_prey_type_changes", &Forager::process_prey_type_changes)
            .def("depletable_detection_probability", &Forager::depletable_detection_probability)
            .def("relative_pursuits_by_position", &Forager::relative_pursuits_by_position)
            .def("relative_pursuits_by_position_single_prey_type", &Forager::relative_pursuits_by_position_single_prey_type)
            .def("proportion_of_detections_within", &Forager::proportion_of_detections_within)
            .def("spatial_detection_proportions", &Forager::spatial_detection_proportions)
            .def("get_prey_type", &Forager::get_prey_type)
            .def("get_favorite_prey_type", &Forager::get_favorite_prey_type)
            .def("get_prey_types", &Forager::get_prey_types)
            .def("get_strategy", &Forager::get_strategy)
            .def("get_parameter", &Forager::get_parameter)
            .def("get_strategy_bounds", &Forager::get_strategy_bounds)
            .def("get_parameter_bounds", &Forager::get_parameter_bounds)
            .def("get_prey_type", &Forager::get_prey_type)
            .def("get_field_of_view", &Forager::get_field_of_view)
            .def("get_fork_length_cm", &Forager::get_fork_length_cm)
            .def("get_max_radius", &Forager::get_max_radius)
            .def("get_focal_velocity", &Forager::get_focal_velocity)
            .def("get_foraging_attempt_rate", &Forager::get_foraging_attempt_rate)
            .def("get_proportion_of_attempts_ingested", &Forager::get_proportion_of_attempts_ingested)
            .def("get_diet_proportion_for_prey_type", &Forager::get_diet_proportion_for_prey_type)
            .def("get_angular_resolution", &Forager::get_angular_resolution)
            .def("print_strategy", &Forager::print_strategy)
            .def("print_parameters", &Forager::print_parameters)
            .def("print_status", &Forager::print_status)
            .def("print_analytics", &Forager::print_analytics)
            .def("analyze_results", &Forager::analyze_results);

    py::enum_<Forager::Strategy>(forager, "Strategy")
            .value("delta_min", Forager::Strategy::s_delta_min)
            .value("sigma_A", Forager::Strategy::s_sigma_A)
            .value("mean_column_velocity", Forager::Strategy::s_mean_column_velocity)
            .value("saccade_time", Forager::Strategy::s_saccade_time)
            .value("discrimination_threshold", Forager::Strategy::s_discrimination_threshold)
            .value("search_image", Forager::Strategy::s_search_image)
            .export_values();

    py::enum_<Forager::Parameter>(forager, "Parameter")
            .value("delta_0", Forager::Parameter::p_delta_0)
            .value("alpha_tau", Forager::Parameter::p_alpha_tau)
            .value("alpha_d", Forager::Parameter::p_alpha_d)
            .value("beta", Forager::Parameter::p_beta)
            .value("A_0", Forager::Parameter::p_A_0)
            .value("t_s_0", Forager::Parameter::p_t_s_0)
            .value("discriminability", Forager::Parameter::p_discriminability)
            .value("flicker_frequency", Forager::Parameter::p_flicker_frequency)
            .value("tau_0", Forager::Parameter::p_tau_0)
            .value("nu", Forager::Parameter::p_nu)
            .export_values();

    /* --------------------- Optimizer class and components ----------------------- */

    py::class_<Optimizer> optimizer(m, "Optimizer");
    optimizer.def(py::init<Forager *, size_t, size_t, bool>())
            .def("optimize_forager", &Optimizer::optimize_forager)
            .def("set_algorithm_options", &Optimizer::set_algorithm_options)
            .def("add_context", &Optimizer::add_context)
            .def("clear_context", &Optimizer::clear_context);


}
