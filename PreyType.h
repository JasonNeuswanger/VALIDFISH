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
    long long int uniqueid; // manually specified by user during input, used for output-sorting and hash keys; type is long long to avoid casting for memoization
    std::string name;
    bool search_image_eligible;
    double length, energy_content, prey_drift_concentration, debris_drift_concentration, crypticity;


    // Attributes that characterize the fish's current interaction with the prey type
    enum SearchImageStatus { search_image_target, search_image_exclusion, no_search_image };
    SearchImageStatus search_image_status;
    double max_visible_distance, rsq, rhosq;
    double prey_pursuit_rate, debris_pursuit_rate, foraging_attempt_rate, proportion_of_attempts_ingested, diet_proportion;

    PreyType(int number, std::string name, double length, double energy_content,
                 double prey_drift_concentration, double debris_drift_concentration, bool search_image_eligible,
                 double crypticity);
    PreyType(PreyType *otherPreyType);

    bool operator< (const PreyType& pt) const; // For sorting

    void compute_details(double fork_length_cm, double theta);

    std::string get_name();
    double get_length();
    double get_diet_proportion();
    double get_max_visible_distance();
    double get_prey_drift_concentration();
    double get_debris_drift_concentration();
    double get_energy_content();

    };


#endif //DRIFTMODELC_PREYTYPE_H
