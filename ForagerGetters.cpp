//
// Created by Jason Neuswanger on 2/27/18.
//

#include "Forager.h"

std::shared_ptr<PreyType> Forager::get_prey_type(std::string name) {
    for (auto & pt : prey_types) {
        if (pt->name == name) {
            return pt;
        }
    }
    printf("WARNING: No prey type exists with name %s, returning null from request to get_prey_type.\n", name.c_str());
    return nullptr;
}

std::vector<std::shared_ptr<PreyType>> Forager::get_prey_types() {
    return prey_types;
}

double Forager::get_fork_length_cm() {
    return fork_length_cm;
}

double Forager::get_max_radius() {
    return max_radius;
}

double Forager::get_field_of_view() {
    return theta;
}

double Forager::get_strategy(Strategy strategy) {
    switch (strategy) {
        case s_delta_min: return delta_min;
        case s_sigma_A: return sigma_A;
        case s_mean_column_velocity: return mean_column_velocity;
        case s_saccade_time: return saccade_time;
        case s_discrimination_threshold: return discrimination_threshold;
        case s_search_image: return search_image;
    }
}

double Forager::get_parameter(Parameter parameter) {
    switch (parameter) {
        case p_delta_0: return delta_0;
        case p_alpha_tau: return alpha_tau;
        case p_alpha_d: return alpha_d;
        case p_beta: return beta;
        case p_A_0: return A_0;
        case p_t_s_0: return t_s_0;
        case p_discriminability: return discriminability;
        case p_flicker_frequency: return flicker_frequency;
        case p_tau_0: return tau_0;
        case p_nu: return nu;
    }
}

double Forager::get_focal_velocity() {
    return focal_velocity;
}

double Forager::get_focal_swimming_cost() {
    return focal_swimming_cost;
}

double Forager::get_foraging_attempt_rate() {
    return foraging_attempt_rate;
}

double Forager::get_proportion_of_attempts_ingested() {
    return proportion_of_attempts_ingested;
}

double Forager::get_diet_proportion_for_prey_type(std::shared_ptr<PreyType> pt) {
    if (pt == nullptr) {
        return 0;   // return proportion of 0 if category isn't in diet at all
    } else {
        return pt->get_diet_proportion();
    }
}

double Forager::get_angular_resolution() {
    return angular_resolution;
}

std::array<double, 2> Forager::get_strategy_bounds(Strategy strategy) {
    return strategy_bounds[strategy];
}

std::array<double, 2> Forager::get_parameter_bounds(Parameter parameter) {
    return parameter_bounds[parameter];
}