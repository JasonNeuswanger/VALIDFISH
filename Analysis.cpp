//
// Created by Jason Neuswanger on 2/8/18.
//

#include "Forager.h"

double Forager::relative_pursuits_by_position_single_prey_type(double x, double y, double z, std::shared_ptr<PreyType> pt) {
    /* This function gives an index of the number of pursuits per unit time made on items detected at each position.
     * However, it's based on probability densities, so it can't really translate into meaningful units unless it's
     * integrated over some volume. */
    const double xsq = gsl_pow_2(x);
    const double ysq = gsl_pow_2(y);
    const double zsq = gsl_pow_2(z);
    if (xsq + ysq + zsq > pt->rsq) return 0;
    const double v = water_velocity(z);
    if (!location_is_profitable(x, y, z, *pt)) { return 0; }
    const double detection_pdf = detection_pdf_at_y(y, x, z, *pt);
    const double t_y = time_at_y(y, x, z, *pt);
    auto dps = discrimination_probabilities(t_y, x, z, *pt);
    const double false_positive_probability = dps.first;
    const double true_hit_probability = dps.second;
    const double answer = detection_pdf * v * (pt->prey_drift_concentration * true_hit_probability + pt->debris_drift_concentration * false_positive_probability);
    return answer;
}

double Forager::relative_pursuits_by_position(double x, double y, double z) {
    /* Adds results from the above function across all prey types */
    if (z >= surface_z || z <= bottom_z) return 0;
    double total = 0;
    for (auto & pt : prey_types) {
        total += relative_pursuits_by_position_single_prey_type(fabs(x), y, z, pt); // Using abs(x) here for symmetry
    }
    return total;
}

double Forager::depleted_prey_concentration_single_prey_type(double x, double y, double z, std::shared_ptr<PreyType> pt) {
    // Estimates the concentration of prey at each given point, taking into account the possibility of having not
    // detected it yet, or having detected it but failed to identify it as prey. Units of items/m3.
    // Todo make this return an intuitive number if outside the reaction volume, i.e. upstream or downstream.
    // Todo adjust for handling time by multiplying by proportion of fish's time not spent handling. (1 / 1 + denominator of NREI)
    const double t_y = time_at_y(y, x, z, *pt);
    const double probability_of_not_being_detected_yet = exp(-mean_value_function(t_y, x, z, *pt));
    const double true_hit_probability = discrimination_probabilities(t_y, x, z, *pt).second;
    return pt->prey_drift_concentration * (probability_of_not_being_detected_yet + (1 - probability_of_not_being_detected_yet) * (1 - true_hit_probability));
}

double Forager::depleted_prey_concentration_total_energy(double x, double y, double z) {
    // Combines depleted concentrations of all prey types to get total energy content remaining.
    double total = 0;
    for (auto & pt : prey_types) {
        total += pt->energy_content * depleted_prey_concentration_single_prey_type(x, y, z, pt);
    }
    return total;
}

std::map<std::string, std::vector<std::map<std::string, double>>> Forager::spatial_detection_proportions(std::shared_ptr<PreyType> pt, std::string which_items, bool verbose) {
    // returns a map with keys 'distance' and 'angle', each of which has as its value a vector of
    // maps, with each entry of the top-level vector containing a map with keys 'min_angle', 'max_angle', 'proportion' or similar for distance
    const double bodylength_m = 0.01 * fork_length_cm;
    std::vector<std::pair<double, double>> distance_bins = {
            {0, 0.5 * bodylength_m},
            {0.5 * bodylength_m, 1.0 * bodylength_m},
            {1.0 * bodylength_m, 2.0 * bodylength_m},
            {2.0 * bodylength_m, 4.0 * bodylength_m},
            {4.0 * bodylength_m, 8.0 * bodylength_m}
    };
    std::vector<std::pair<double, double>> angle_bins = {
            {0, 0.25 * M_PI},
            {0.25 * M_PI, 0.5 * M_PI},
            {0.5 * M_PI, 0.75 * M_PI},
            {0.75 * M_PI, M_PI},
            {M_PI, 1.25 * M_PI},
            {1.25 * M_PI, 1.5 * M_PI}
    };
    std::vector<std::shared_ptr<PreyType>> types;
    if (pt == nullptr) {
        types = prey_types;                              // Run through all types if type was given as "All"
    } else {
        types.push_back(pt);
    }
    auto main_lambda = [this, which_items](std::shared_ptr<PreyType> ipt, double min_distance, double max_distance, double min_angle, double max_angle) {  // integrates same function as above, but over the whole search volume
        if (ipt->prey_drift_concentration == 0 && ipt->debris_drift_concentration == 0) { return 0.0; }
        double inner_theta, inner_phi;  // placeholders for the inner integrands
        double max_possible_distance = ipt->get_max_visible_distance();
        double truncated_min_distance = (min_distance > max_possible_distance) ? max_possible_distance : min_distance;
        double truncated_max_distance = (max_distance > max_possible_distance) ? max_possible_distance : max_distance;
        if (truncated_min_distance == truncated_max_distance) { return 0.0; }
        auto integrand = [this, ipt, which_items](double rho, double *theta, double *phi)->double{
            cartesian_3D_coords coords = cartesian_from_spherical(rho, *theta, *phi);
            if (coords.z > surface_z || coords.z < bottom_z) { return 0; }
            if (!location_is_profitable(coords.x, coords.y, coords.z, *ipt)) { return 0; }
            double v = water_velocity(coords.z);
            double prob = detection_pdf_at_y(coords.y, coords.x, coords.z, *ipt) * v;
            double t = (sqrt(ipt->rsq - gsl_pow_2(coords.x) - gsl_pow_2(coords.z)) - coords.y) / v;
            auto dps = discrimination_probabilities(t, coords.x, coords.z, *ipt);
            double false_positive_probability = dps.first;
            double true_hit_probability = dps.second;
            if (which_items == "Prey") {
                prob *= ipt->prey_drift_concentration * true_hit_probability;
            } else if (which_items == "Debris") {
                prob *= ipt->debris_drift_concentration * false_positive_probability;
            } else if (which_items == "All") {
                prob *= (ipt->prey_drift_concentration * true_hit_probability + ipt->debris_drift_concentration * false_positive_probability);
            }
            if (isnan(prob)) {
                printf("******** NaN returned for totals at v=%.8f, (x,y,z)=(%.5f, %.5f, %.5f), for prey type %s items %s.\n", v, coords.x, coords.y, coords.z, ipt->name.c_str(), which_items.c_str());
            }
            return prob;
        };
        gsl_function_pp_3d<decltype(integrand)> Fintegrand(integrand, &inner_theta, &inner_phi);
        gsl_function *F = static_cast<gsl_function*>(&Fintegrand);
        return integrate_over_volume(F, truncated_min_distance, truncated_max_distance, min_angle, max_angle);
    };
    if (verbose) {
        printf("Computing distance/angle proportions...\n");
    }
    double total = 0;
    for (auto &jpt : types) {
        total += main_lambda(jpt, NAN, NAN, NAN, NAN);
    }
    std::vector<std::map<std::string, double>> distance_results, angle_results;
    for (auto distance_pair : distance_bins) {
        double min_distance = distance_pair.first;
        double max_distance = distance_pair.second;
        double bin_total = 0;
        for (auto & jpt : types) {
            bin_total += main_lambda(jpt, min_distance, max_distance, NAN, NAN);
        }
        std::map<std::string, double> bin_result = {
                {"min_distance", min_distance},
                {"max_distance", max_distance},
                {"proportion", bin_total / total}
        };
        distance_results.push_back(bin_result);
        if (verbose) {
            printf("For distance bin %.3f to %.3f m,    predicted proportion is %.3f.\n", min_distance, max_distance, bin_total / total);
        }
    }
    for (auto angle_pair : angle_bins) {
        double min_angle = angle_pair.first;
        double max_angle = angle_pair.second;
        double bin_total = 0;
        for (auto & jpt : types) {
            bin_total += main_lambda(jpt, NAN, NAN, min_angle, max_angle);
        }
        std::map<std::string, double> bin_result = {
                {"min_angle", min_angle},
                {"max_angle", max_angle},
                {"proportion", bin_total / total}
        };
        angle_results.push_back(bin_result);
        if (verbose) {
            printf("For angle bin %.3f to %.3f radians, predicted proportion is %.3f.\n", min_angle, max_angle, bin_total / total);
        }
    }
    std::map<std::string, std::vector<std::map<std::string, double>>> combined_results = {
            {"distance", distance_results},
            {"angle", angle_results}
    };
    return combined_results;
}

std::shared_ptr<PreyType> Forager::get_favorite_prey_type() {
    // Returns the favorite prey type by quantity of prey pursued
    std::shared_ptr<PreyType> favorite;
    double favorite_pursuit_rate = -1;
    for (auto & pt : prey_types) {
        double rate = pursuit_rate("prey", pt);
        if (rate > favorite_pursuit_rate) {
            favorite = pt;
            favorite_pursuit_rate = rate;
        }
    }
    return favorite;
}

double Forager::pursuit_rate(std::string which_rate, std::shared_ptr<PreyType> pt) {
    /* For which_rate, pass either "prey" or "debris." For prey type, pass a specific prey type object or nullptr
     * to sum over all prey types. */
    double x, y; // temporary holder for x values (y is required but unused) during the integrations
    assert(which_rate == "prey" || which_rate == "debris");
    if (pt != nullptr && ((which_rate == "prey" && pt->prey_drift_concentration == 0) || (which_rate == "debris" && pt->debris_drift_concentration == 0))) {
        return 0;
    }
    std::vector<std::shared_ptr<PreyType>> types;
    if (pt == nullptr) {
        types = prey_types;                              // Run through all types if type was given as NULL
    } else {
        types.push_back(pt);
    }
    auto numerator_integrand = [this, types, which_rate](double z, double *y, double *x) -> double {
        double sum = 0;
        const double v = water_velocity(z);
        bool prey_is_resolvable;
        double pd, pf, ph, dp, dd;
        for (auto & ipt : types) {
            prey_is_resolvable = (gsl_pow_2(*x) + gsl_pow_2(z) <= ipt->rsq);
            if (prey_is_resolvable && (ipt->prey_drift_concentration > 0 || ipt->debris_drift_concentration > 0)) {
                pd = detection_probability(*x, z, *ipt);
                auto mdps = expected_discrimination_probabilities(*x, z, *ipt, pd);
                pf = mdps.first;     // false positive probability
                ph = mdps.second;    // true hit probability
                dp = ipt->prey_drift_concentration;
                dd = ipt->debris_drift_concentration;
                if (which_rate == "prey") {
                    sum += (pd * v * (ph * dp));
                } else if (which_rate == "debris") {
                    sum += (pd * v * (pf * dd));
                }
            }
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(numerator_integrand)> Fn(numerator_integrand, &y, &x);
    gsl_function *Fnumerator = static_cast<gsl_function *>(&Fn);
    double numerator = integrate_over_xz_plane(Fnumerator, false);
    auto denominator_integrand = [this](double z, double *y, double *x) -> double {
        double sum = 0;
        const double v = water_velocity(z);
        bool prey_is_resolvable;
        double pd, pf, ph, h, dp, dd;
        for (auto &ipt : prey_types) {
            prey_is_resolvable = (gsl_pow_2(*x) + gsl_pow_2(z) <= ipt->rsq);
            if (prey_is_resolvable && (ipt->prey_drift_concentration > 0 || ipt->debris_drift_concentration > 0)) {
                pd = detection_probability(*x, z, *ipt);
                auto mdps = expected_discrimination_probabilities(*x, z, *ipt, pd);
                pf = mdps.first;     // false positive probability
                ph = mdps.second;    // true hit probability
                h = expected_maneuver_cost(*x, z, *ipt, false, pd);
                dp = ipt->prey_drift_concentration;
                dd = ipt->debris_drift_concentration;
                sum += pd * v * h * (ph * dp + pf * dd);
            }
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(denominator_integrand)> Fd(denominator_integrand, &y, &x);
    gsl_function *Fdenominator = static_cast<gsl_function *>(&Fd);
    double denominator = integrate_over_xz_plane(Fdenominator, false);
    double rate = numerator / (1 + denominator);
    if (isnan(rate) && pt != nullptr) {
        printf("ERROR: Somehow got NaN pursuit rate for prey type %s with numerator %.4f and denominator %.4f.\n", pt->name.c_str(), numerator, denominator);
    } else if (isnan(rate)) {
        printf("ERROR: Somehow got NaN pursuit rate for prey type null (all types) with numerator %.4f and denominator %.4f.\n", numerator, denominator);
    }
    assert(!isnan(rate));
    return rate;
}

void Forager::analyze_results() {
    /* Master command to compute all the metrics relevant for comparing the model to the data, following an optimization. */
    prey_pursuit_rate = pursuit_rate("prey", nullptr);
    debris_pursuit_rate = pursuit_rate("debris", nullptr);
    foraging_attempt_rate = prey_pursuit_rate + debris_pursuit_rate;
    proportion_of_attempts_ingested = prey_pursuit_rate / foraging_attempt_rate;
    for (auto & pt : prey_types) {
        pt->prey_pursuit_rate = pursuit_rate("prey", pt);
        pt->debris_pursuit_rate = pursuit_rate("debris", pt);
        pt->foraging_attempt_rate = pt->prey_pursuit_rate + pt->debris_pursuit_rate;
        pt->proportion_of_attempts_ingested = pt->prey_pursuit_rate / pt->foraging_attempt_rate;
        pt->diet_proportion = pt->prey_pursuit_rate / prey_pursuit_rate;
    }
}

