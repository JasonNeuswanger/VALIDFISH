//
// Created by Jason Neuswanger on 2/1/18.
//

#ifndef DRIFTMODELC_PREYCATEGORY_H
#define DRIFTMODELC_PREYCATEGORY_H

#include <string>
#include <cmath>

class PreyCategory {

private:

    //static short category_count;    // automaticaly incremented unique index for sorting/hashing prey categories

    double feature_alpha, attention_allocated;

public:

    int uniqueid;                    // instance variable set on initialization based on category_count

    std::string name;

    double length, energy_content, prey_drift_density, debris_drift_density, feature_size, crypticity;

    double false_positive_probability, true_hit_probability;

    double prey_pursuit_rate, debris_pursuit_rate, foraging_attempt_rate, proportion_of_attempts_ingested, diet_proportion;

    PreyCategory(int number,
                   std::string _name,
                   double _length,                  // length in m
                   double _energy_content,          // energy content in J, or -1 if unknown
                   double _prey_drift_density,      // prey drift density in items/m3
                   double _debris_drift_density,    // debris drift density in items/m3
                   double _feature_size,            // size as a proportion of feature space
                   double _crypticity);
    PreyCategory(PreyCategory *otherPreyCategory);

    bool operator< (const PreyCategory& pc) const; // For sorting

    double get_attention_allocated();
    void set_attention_allocated(double attention);
    double get_feature_alpha();
    double get_diet_proportion();
    double max_visible_distance(double fish_fork_length_cm);

};


#endif //DRIFTMODELC_PREYCATEGORY_H
