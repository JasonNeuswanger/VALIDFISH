//
// Created by Jason Neuswanger on 2/1/18.
//

#include "PreyCategory.h"

//short PreyCategory::category_count = 0;

PreyCategory::PreyCategory(int number,                    // Unique number for sorting
                           std::string _name,               // Main label
                           double _length,                  // length in m
                           double _energy_content,          // energy content in J, or -1 if unknown
                           double _prey_drift_density,      // prey drift density in items/m3
                           double _debris_drift_density,    // debris drift density in items/m3
                           double _feature_size,            // size as a proportion of feature space
                           double _crypticity) {            // unitless crypticity
    //uniqueid = ++category_count;
    uniqueid = number;
    name = _name;
    length = _length;
    prey_drift_density = _prey_drift_density;
    debris_drift_density = _debris_drift_density;
    feature_size = _feature_size;
    crypticity = _crypticity;
    if (_energy_content == -1) {
        double dry_mass = exp(-5.021 + 2.88 * log(1000 * length));  // dry mass in mg (Smock 1980)
        energy_content = 28.3 * dry_mass;                           //  energy in J (Cummins & Wuycheck 1971)
    } else {
        energy_content = _energy_content;
    }
    set_attention_allocated(0);
}

PreyCategory::PreyCategory(PreyCategory *otherPreyCategory) {
    /* Creates a deep copy of another prey category, for use deep copying fish */
    //uniqueid = ++category_count;
    uniqueid = otherPreyCategory->uniqueid;
    name = otherPreyCategory->name;
    length = otherPreyCategory->length;
    prey_drift_density = otherPreyCategory->prey_drift_density;
    debris_drift_density = otherPreyCategory->debris_drift_density;
    feature_size = otherPreyCategory->feature_size;
    crypticity = otherPreyCategory->crypticity;
    energy_content = otherPreyCategory->energy_content;
    set_attention_allocated(otherPreyCategory->attention_allocated);
}

bool PreyCategory::operator< (const PreyCategory& pc) const
{
    return (uniqueid < pc.uniqueid);
}

double PreyCategory::get_attention_allocated() {
    return attention_allocated;
}

void PreyCategory::set_attention_allocated(double attention) {
    attention_allocated = attention;
    feature_alpha = attention_allocated / feature_size;
}

double PreyCategory::get_feature_alpha() {
    return feature_alpha;
}

double PreyCategory::get_diet_proportion() {
    return diet_proportion;
}

double PreyCategory::max_visible_distance(double fish_fork_length_cm) {
    /* From Hughes & Dill 2000, but with units converted so prey length and the returned answer are both in m */
    return 120. * length * (1. - exp(-0.2 * fish_fork_length_cm));
}