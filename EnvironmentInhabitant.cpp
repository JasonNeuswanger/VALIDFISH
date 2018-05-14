//
// Created by Jason Neuswanger on 1/31/18.
//

#include "Forager.h"

void Forager::compute_focal_velocity() {
    focal_velocity = water_velocity(0);
}

//double Forager::water_velocity(double z) {
//    /* I tested water velocity below with a cache, but the cache hits took 30-50 % longer than just recalculating. */
//    double k = bed_roughness * 1000;       // bed roughness height in mm
//    double R = depth;                      // hydraulic delta_min in m, approximated as depth, as per Hayes et al 2007
//    double H = 0.001 + surface_z - z;      //  distance of z position below the surface, in m; adding 1 mm to avoid infinite velocity at surface
//    double vstar = mean_column_velocity / (5.75 * log10(12.27 * R / k));  // # Stream Hydrology: An Introduction for Ecologists eqn 6.50
//    return 5.75 * log10(30 * H / k) * vstar;  // Hayes et al 2007 eqn 1
//}

double Forager::water_velocity(double z) {
    return 1.4 * mean_column_velocity * pow((z-bottom_z)/depth, 0.52);  // todo get actual constants for this
}

void Forager::build_sample_prey_types() {
//    add_prey_type(1, "1 mm class",   0.001,  -1, 1, 30,    50000, false);
//    add_prey_type(2, "2 mm class",   0.002,  -1, 1, 9,     3000,  false);
//    add_prey_type(3, "3-4 mm class", 0.0035, -1, 1, 4,     700,   false);
//    add_prey_type(4, "5-7 mm class", 0.006,  -1, 1, 0.15,  125,   false);
//    add_prey_type(5, "8+ mm class",  0.01,   -1, 1, 0.015, 30,    false);

    // The problem is that the cache is frequently missing when I use all these prey types. However, when I
    // use the sample prey types with all the same fish, the caches often hit.
    // Ahh, units on prey length are one of the problems.

    add_prey_type(5, "5 mm size class", 0.0050, 9.6548, 1, 0.1948, 2.3182, false);
    add_prey_type(10, "Cinygmula adult", 0.0000, 0.0000, 1, 0.0000, 0.0000, false);
    add_prey_type(8, "10+ mm size class", 0.0100, 4.2332, 1, 0.0600, 0.1818, false);
    add_prey_type(9, "Drunella adult", 0.0000, 0.0000, 1, 0.0000, 0.0000, false);
    add_prey_type(7, "8-9 mm size class", 0.0090, 3.0058, 1, 0.0600, 0.1475, false);
    add_prey_type(6, "6-7 mm size class", 0.0063, 15.3199, 1, 0.0899, 0.5339, false);
    add_prey_type(4, "4 mm size class", 0.0040, 2.3011, 1, 0.8918, 5.2908, false);
    add_prey_type(2, "2 mm size class", 0.0020, 0.9435, 1, 8.1762, 101.5793, false);
    add_prey_type(1, "1 mm size class", 0.0010, 0.3865, 1, 7.7640, 2240.3654, false);
    add_prey_type(3, "3 mm size class", 0.0030, 2.8867, 1, 2.1284, 18.0124, false);
    process_prey_type_changes();
}

void Forager::add_prey_type(int number, std::string name, double mean_prey_length, double mean_prey_energy,
                            double crypticity, double prey_drift_concentration, double debris_drift_concentration,
                            bool search_image_eligible) {
    auto pt = std::make_shared<PreyType>(number, name, mean_prey_length, mean_prey_energy, prey_drift_concentration,
                                         debris_drift_concentration, search_image_eligible, crypticity);
    if ((pt->length < min_prey_length_from_gill_rakers) || (pt->length > max_prey_length_from_mouth_gape)) {
        // If a prey type is too small or too large to physically eat, treat all prey items within that type as debris.
        pt->debris_drift_concentration += pt->prey_drift_concentration;
        pt->prey_drift_concentration = 0;
    }
    prey_types.push_back(pt);
}

void Forager::process_prey_type_changes() {
    std::sort(prey_types.begin(), prey_types.end());
    // How mapping search image value (a double between -1 and 1) onto prey types works:
    // This is a bit messy because it's a categorical variable but optimization with the Grey Wolf algorithm will
    // work best if it's represented as a single continuous variable. Therefore, we map values from -1 to 0 onto
    // a "no search image" strategy. Values from 0 to 1 are divided at even intervals among the prey types eligible
    // for search images. This gives the algorithm the best chance of considering all the strategies.
    std::vector<std::shared_ptr<PreyType>> search_image_eligible_prey_types;
    for (auto & pt : prey_types) {
        if (pt->search_image_eligible) {
            search_image_eligible_prey_types.push_back(pt);
        }
    }
    size_t n_eligible_types = search_image_eligible_prey_types.size();
    if (search_image > 0 && n_eligible_types > 0) {
        for (auto & pt : prey_types) {
            pt->search_image_status = PreyType::SearchImageStatus::search_image_exclusion;   // exclude all types
        }
        double threshold_increment = 1.0 / n_eligible_types;                                // then include the right one
        for (int i=0; i<n_eligible_types; i++) {
            if (i*threshold_increment < search_image && search_image < (i+1)*threshold_increment) {
                search_image_eligible_prey_types[i]->search_image_status = PreyType::SearchImageStatus::search_image_target;
            }
        }
    } else {
        for (auto & pt : prey_types) {
            pt->search_image_status = PreyType::SearchImageStatus::no_search_image;
        }
    }
    process_parameter_updates();
}


size_t Forager::num_prey_types() {
    return prey_types.size();
}

size_t Forager::num_search_image_eligible_prey_types() {
    size_t count = 0;
    for (auto & pt : prey_types) {
        if (pt->search_image_eligible) {
            count++;
        }
    }
    return count;
}
