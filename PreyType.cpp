//
// Created by Jason Neuswanger on 2/1/18.
//

#include "Forager.h"
#include "PreyType.h"

PreyType::PreyType(int number, std::string name, double length, double energy_content,
                   double prey_drift_concentration, double debris_drift_concentration, bool search_image_eligible,
                   double crypticity) {                 // unitless crypticity
    this->uniqueid = number;
    this->name = std::move(name);
    this->length = length;
    this->real_prey_drift_concentration = prey_drift_concentration;     // Holders to retain the real values so the main variables
    this->real_debris_drift_concentration = debris_drift_concentration; // can be manipulated to plot response to changes, etc.
    this->prey_drift_concentration = prey_drift_concentration;
    this->debris_drift_concentration = debris_drift_concentration;
    this->search_image_eligible = search_image_eligible;
    this->crypticity = crypticity;
    if (energy_content == -1) {
        double dry_mass = exp(-5.021 + 2.88 * log(1000 * length));  // dry mass in mg (Smock 1980)
        this->energy_content = 28.3 * dry_mass;                     //  energy in J (Cummins & Wuycheck 1971)
    } else {
        this->energy_content = energy_content;
    }
    search_image_status = no_search_image;
}

PreyType::PreyType(PreyType *otherPreyType) {
    /* Creates a deep copy of another prey type, for use deep copying fish */
    uniqueid = otherPreyType->uniqueid;
    name = otherPreyType->name;
    length = otherPreyType->length;
    prey_drift_concentration = otherPreyType->prey_drift_concentration;
    debris_drift_concentration = otherPreyType->debris_drift_concentration;
    search_image_eligible = otherPreyType->search_image_eligible;
    search_image_status = otherPreyType->search_image_status;
    crypticity = otherPreyType->crypticity;
    energy_content = otherPreyType->energy_content;
}

bool PreyType::operator< (const PreyType& pt) const
{
    return (uniqueid < pt.uniqueid);
}

std::string PreyType::get_name() {
    return name;
}

double PreyType::get_length() {
    return length;
}

void PreyType::set_prey_concentration_multiplier(double value) {
    prey_concentration_multiplier = value;
    prey_drift_concentration = real_prey_drift_concentration * prey_concentration_multiplier;
}

void PreyType::set_debris_concentration_multiplier(double value) {
    debris_concentration_multiplier = value;
    debris_concentration_multiplier = real_debris_drift_concentration * debris_concentration_multiplier;
}

double PreyType::get_prey_drift_concentration() {
    return prey_drift_concentration;
}

double PreyType::get_debris_drift_concentration() {
    return debris_drift_concentration;
}

double PreyType::get_energy_content() {
    return energy_content;
}

double PreyType::get_diet_proportion() {    // Only valid after calling forager->analyze_results()
    return diet_proportion;
}

double PreyType::get_prey_pursuit_rate() {
    return prey_pursuit_rate;
}

double PreyType::get_debris_pursuit_rate() {
    return debris_pursuit_rate;
}


double PreyType::get_max_visible_distance() {
    return max_visible_distance;
}

void PreyType::compute_details(double fork_length_cm, double theta) {
    // Computes details of the prey type that depend on passed attributes of the fish.
    max_visible_distance = 120. * length * (1. - exp(-0.2 * fork_length_cm));  // Hughes & Dill 1990, but with units converted so prey length and the returned answer are both in m
    rsq = gsl_pow_2(max_visible_distance); // Shorthand to cut a few computations elsewhere
    rhosq = (theta < M_PI) ? gsl_pow_2(max_visible_distance * sin(theta / 2)) : rsq; // rho = max prey_radius in lateral direction
//    if (isnan(rhosq)) {
//        printf("WARNING: rhosq is NaN for prey type %s for fish with fork length = %.2f cm and theta = %.2f radians.\n", name.c_str(), fork_length_cm, theta);
//    } else if (rhosq == 0) {
//        printf("WARNING: rhosq is 0 for prey type %s for fish with fork length = %.2f cm and theta = %.2f radians.\n", name.c_str(), fork_length_cm, theta);
//    } else {
//        printf("Setting rhosq = %.8f for prey type %s for fish with fork length = %.2f cm and theta = %.2f radians.\n", rhosq, name.c_str(), fork_length_cm, theta);
//    }
//    assert(!isnan(rhosq));
//    assert(rhosq > 0);
}