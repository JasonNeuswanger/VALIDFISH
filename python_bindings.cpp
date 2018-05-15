//
// Created by Jason Neuswanger on 2/21/18.
//

#include "lib/pybind11/include/pybind11/pybind11.h"
#include "lib/pybind11/include/pybind11/stl.h"
#include "Forager.h"
#include "Optimizer.h"

namespace py = pybind11;

PYBIND11_MODULE(pyvalidfish, m) {

    m.doc() = "Visual Attention-Limited Individual Drift Forager In Simulated Habitat (VALIDFISH) model in C++.";

    /* --------------------- Prey type class and components ----------------------- */

    py::class_<PreyType, std::shared_ptr<PreyType>> preytype(m, "PreyType");

    preytype.def(py::init<int, std::string, double, double, double, double, bool, double>())
            .def("get_diet_proportion", &PreyType::get_diet_proportion)
            .def("get_max_visible_distance", &PreyType::get_max_visible_distance)
            .def("get_name", &PreyType::get_name)
            .def("get_length", &PreyType::get_length)
            .def("get_prey_drift_concentration", &PreyType::get_prey_drift_concentration)
            .def("get_debris_drift_concentration", &PreyType::get_debris_drift_concentration)
            .def("get_energy_content", &PreyType::get_energy_content);

    /* --------------------- Forager class and components ----------------------- */

    py::class_<Forager, std::shared_ptr<Forager>> forager(m, "Forager");

    forager.def(py::init<double, double, double, double, double, double, double, double, double, double, double, double,
                    double, double, unsigned, double, double, double, double , double, double, double, double, std::string *>())
            .def("NREI", &Forager::NREI)
            .def("GREI", &Forager::GREI)
            .def("maneuver_cost_rate", &Forager::maneuver_cost_rate)
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
            .def("water_velocity", &Forager::water_velocity)
            .def("relative_pursuits_by_position", &Forager::relative_pursuits_by_position)
            .def("relative_pursuits_by_position_single_prey_type", &Forager::relative_pursuits_by_position_single_prey_type)
            .def("depleted_prey_concentration_total_energy", &Forager::depleted_prey_concentration_total_energy)
            .def("depleted_prey_concentration_single_prey_type", &Forager::depleted_prey_concentration_single_prey_type)
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
            .def("get_focal_swimming_cost", &Forager::get_focal_swimming_cost)
            .def("get_foraging_attempt_rate", &Forager::get_foraging_attempt_rate)
            .def("get_proportion_of_attempts_ingested", &Forager::get_proportion_of_attempts_ingested)
            .def("get_diet_proportion_for_prey_type", &Forager::get_diet_proportion_for_prey_type)
            .def("get_angular_resolution", &Forager::get_angular_resolution)
            .def("print_strategy", &Forager::print_strategy)
            .def("print_parameters", &Forager::print_parameters)
            .def("print_status", &Forager::print_status)
            .def("time_NREIs", &Forager::time_NREIs)
            .def("print_analytics", &Forager::print_analytics)
            .def("analyze_results", &Forager::analyze_results);

    py::enum_<Forager::Strategy>(forager, "Strategy")
            .value("sigma_A", Forager::Strategy::s_sigma_A)
            .value("mean_column_velocity", Forager::Strategy::s_mean_column_velocity)
            .value("inspection_time", Forager::Strategy::s_inspection_time)
            .value("discrimination_threshold", Forager::Strategy::s_discrimination_threshold)
            .value("search_image", Forager::Strategy::s_search_image)
            .export_values();

    py::enum_<Forager::Parameter>(forager, "Parameter")
            .value("delta_0", Forager::Parameter::p_delta_0)
            .value("alpha_tau", Forager::Parameter::p_alpha_tau)
            .value("alpha_d", Forager::Parameter::p_alpha_d)
            .value("beta", Forager::Parameter::p_beta)
            .value("A_0", Forager::Parameter::p_A_0)
            .value("discriminability", Forager::Parameter::p_discriminability)
            .value("flicker_frequency", Forager::Parameter::p_flicker_frequency)
            .value("tau_0", Forager::Parameter::p_tau_0)
            .value("nu_0", Forager::Parameter::p_nu_0)
            .value("delta_p", Forager::Parameter::p_delta_p)
            .value("omega_p", Forager::Parameter::p_omega_p)
            .value("ti_p", Forager::Parameter::p_ti_p)
            .value("sigma_p_0", Forager::Parameter::p_sigma_p_0)
            .export_values();

    /* --------------------- Optimizer class and components ----------------------- */

    py::class_<Optimizer, std::shared_ptr<Optimizer>> optimizer(m, "Optimizer");
    optimizer.def(py::init<std::shared_ptr<Forager>, size_t, size_t, bool>())
            .def("optimize_forager", &Optimizer::optimize_forager)
            .def("add_context", &Optimizer::add_context)
            .def("clear_context", &Optimizer::clear_context);


}
