//
// Created by Jason Neuswanger on 1/31/18.
//

#include "Forager.h"

void Forager::compute_focal_velocity() {
    focal_velocity = water_velocity(0);
}

double Forager::water_velocity(double z) {
    /* I tested water velocity below with a cache, but the cache hits took 30-50 % longer than just recalculating. */
    double k = bed_roughness * 1000;       // bed roughness height in mm
    double R = depth;                      // hydraulic radius in m, approximated as depth, as per Hayes et al 2007
    double H = 0.001 + surface_z - z;      //  distance of z position below the surface, in m; adding 1 mm to avoid infinite velocity at surface
    double vstar = mean_column_velocity / (5.75 * log10(12.27 * R / k));  // # Stream Hydrology: An Introduction for Ecologists eqn 6.50
    return 5.75 * log10(30 * H / k) * vstar;  // Hayes et al 2007 eqn 1
}

void Forager::normalize_feature_sizes() {
    double total_feature_size = 0;
    for (auto & pc : prey_categories) {
        total_feature_size += pc.feature_size;
    }
    for (auto & pc : prey_categories) {
        pc.feature_size = pc.feature_size / total_feature_size;
    }
}

void Forager::evenly_distribute_attention() {
    for (auto & pc : prey_categories) {
        pc.set_attention_allocated(pc.feature_size);
    }
}

void Forager::build_sample_prey_categories() {
    add_prey_category(1, "1 mm class", 0.001, -1, 1.0, 30, 50000, 0.2);
    add_prey_category(2, "2 mm class", 0.002, -1, 1.0, 9, 3000, 0.2);
    add_prey_category(3, "3-4 mm class", 0.0035, -1, 1.0, 4, 700, 0.2);
    add_prey_category(4, "5-7 mm class", 0.006, -1, 1.0, 0.15, 125, 0.2);
    add_prey_category(5, "8+ mm class", 0.01, -1, 1.0, 0.015, 30, 0.2);
    process_prey_category_changes();
}

void Forager::add_prey_category(int number, std::string name, double mean_prey_length, double mean_prey_energy,
                                double crypticity_multiplier, double prey_drift_density, double debris_drift_density, double feature_size) {
    prey_categories.emplace_back(PreyCategory(number, name, mean_prey_length, mean_prey_energy, prey_drift_density,
                                              debris_drift_density, feature_size, base_crypticity * crypticity_multiplier));
}

void Forager::process_prey_category_changes() {
    normalize_feature_sizes();
    evenly_distribute_attention();
    std::sort(prey_categories.begin(), prey_categories.end());
    process_parameter_updates();
}

size_t Forager::num_prey_categories() {
    return prey_categories.size();
}
