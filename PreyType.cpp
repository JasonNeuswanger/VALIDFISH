//
// Created by Jason Neuswanger on 2/1/18.
//

#include "Forager.h"
#include "PreyType.h"

PreyType::PreyType(int number, std::string name, double length, double energy_content,
                   double prey_drift_concentration, double debris_drift_concentration, bool search_image_eligible,
                   double crypticity) {                 // unitless crypticity
    this->uniqueid = number;
    this->name = name;
    this->length = length;
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

bool PreyType::operator< (const PreyType& pc) const
{
    return (uniqueid < pc.uniqueid);
}

std::string PreyType::get_name() {
    return name;
}

double PreyType::get_length() {
    return length;
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

double PreyType::get_diet_proportion() {
    return diet_proportion;
}

double PreyType::get_perceptual_sigma() {
    return perceptual_sigma;
}

double PreyType::get_false_positive_probability() {
    return false_positive_probability;
}

double PreyType::get_true_hit_probability() {
    return true_hit_probability;
}

double PreyType::get_max_visible_distance() {
    return max_visible_distance;
}

double PreyType::get_max_attended_distance() {
    return max_attended_distance;
}

void PreyType::compute_details(double fork_length_cm, double saccade_time, double t_s_0, double discrimination_threshold, double discriminability, double alpha_d, double delta_min) {
    max_visible_distance = 120. * length * (1. - exp(-0.2 * fork_length_cm));  // Hughes & Dill 1990, but with units converted so prey length and the returned answer are both in m
    max_attended_distance = fmin(max_visible_distance, length / (M_PI * tan(delta_min / 2)));
    double search_image_effect = 1;
    switch (search_image_status) {
        case no_search_image:
            break; // stick with the default value of 1
        case search_image_exclusion:
            search_image_effect = 10000;
            break;
        case search_image_target:
            search_image_effect = alpha_d;
            break;
    }
    perceptual_sigma = sqrt(1 + t_s_0 / (search_image_effect * saccade_time));

    false_positive_probability = 1 - gsl_cdf_gaussian_P(discrimination_threshold / perceptual_sigma, 1);
    true_hit_probability = 1 - gsl_cdf_gaussian_P((discrimination_threshold - discriminability) / perceptual_sigma, 1);

    probability_a_pursued_item_is_prey = prey_drift_concentration * true_hit_probability /
                                                (prey_drift_concentration * true_hit_probability + debris_drift_concentration * false_positive_probability);
    expected_energy_gain = energy_content * probability_a_pursued_item_is_prey;
}