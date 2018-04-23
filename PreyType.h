//
// Created by Jason Neuswanger on 2/1/18.
//

#ifndef DRIFTMODELC_PREYTYPE_H
#define DRIFTMODELC_PREYTYPE_H

#include <string>
#include <cmath>

class PreyType {

public:

    // Attributes that characterize the prey type
    int uniqueid; // manually specified by user during input, used for output-sorting purposes
    std::string name;
    bool search_image_eligible;
    double length, energy_content, prey_drift_concentration, debris_drift_concentration, crypticity;


    // Attributes that characterize the fish's current interaction with the prey type
    enum SearchImageStatus { search_image_target, search_image_exclusion, no_search_image };
    SearchImageStatus search_image_status;
    double perceptual_sigma, false_positive_probability, true_hit_probability, max_visible_distance, max_attended_distance;
    double probability_a_pursued_item_is_prey, expected_energy_gain; // expected energy gain factors in the probability the item will be prey

    double prey_pursuit_rate, debris_pursuit_rate, foraging_attempt_rate, proportion_of_attempts_ingested, diet_proportion;

    PreyType(int number, std::string name, double length, double energy_content,
                 double prey_drift_concentration, double debris_drift_concentration, bool search_image_eligible,
                 double crypticity);
    PreyType(PreyType *otherPreyType);

    bool operator< (const PreyType& pc) const; // For sorting

    void compute_details(double fork_length_cm, double saccade_time, double t_s_0, double discrimination_threshold, double discriminability, double alpha_d, double delta_min);

    std::string get_name();
    double get_length();
    double get_diet_proportion();
    double get_perceptual_sigma();
    double get_false_positive_probability();
    double get_true_hit_probability();
    double get_max_visible_distance();
    double get_max_attended_distance();

};


#endif //DRIFTMODELC_PREYTYPE_H
