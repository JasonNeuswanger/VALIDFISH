//
// Created by Jason Neuswanger on 2/27/18.
//

#include "Forager.h"

PreyCategory *Forager::get_prey_category(std::string *name) {
    for (auto & pc : prey_categories) {
        if (pc.name == *name) {
            return &pc;
        }
    }
    return nullptr;
//    fprintf(stderr, "No prey category exists with name %s. Returning first prey category.\n", (*name).c_str());
//    return &prey_categories.at(0);
}

double Forager::get_fork_length_cm() {
    return fork_length_cm;
}

double Forager::get_radius() {
    return radius;
}

double Forager::get_theta() {
    return theta;
}

double Forager::get_focal_velocity() {
    return focal_velocity;
}

double Forager::get_foraging_attempt_rate() {
    return foraging_attempt_rate;
}

double Forager::get_proportion_of_attempts_ingested() {
    return proportion_of_attempts_ingested;
}

double Forager::get_diet_proportion_for_prey_category(std::string *category_name) {
    PreyCategory *pc = get_prey_category(category_name);
    if (pc == nullptr) {
        return 0;   // return proportion of 0 if category isn't in diet at all
    } else {
        return pc->get_diet_proportion();
    }
}

double Forager::get_angular_resolution() {
    return angular_resolution;
}