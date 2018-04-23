//
// Created by Jason Neuswanger on 2/8/18.
//

#include "Forager.h"

double Forager::depletable_detection_probability(double x, double y, double z, PreyType *pt) {
    /* This function gives the probability of getting to a position undetected times the detection PDF at that
     * point. */
    if (z > surface_z || z < bottom_z) { return 0; }
    double prey_radius = pt->get_max_attended_distance();
    const double xsq = gsl_pow_2(x);
    const double ysq = gsl_pow_2(y);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(prey_radius);
    if (xsq + ysq + zsq > rsq) { return 0; }
    const double y0 = sqrt(rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    if (y < yT || y > y0) { return 0; }
    const double v = water_velocity(z);
    const double t_y = (y0 - y) / v; // Time at which the particle reaches position y
    auto integrand = [this, x, z, pt](double t)->double{
        double tauval = tau(t, x, z, pt);
        return exp(-t / tauval) / tauval;
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
        gsl_integration_qags(F, 0, t_y, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
        gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, 0, t_y, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    const double probability_of_no_previous_detection = 1 - result;
    const double tau_y = tau(t_y, x, z, pt);
    const double probability_of_detection_at_y = exp(-t_y / tau_y) / tau_y; // probability of detecting in one second (detection PDF) at position y
    const double answer = probability_of_no_previous_detection * probability_of_detection_at_y;
    return answer;
}

double Forager::relative_pursuits_by_position_single_prey_type(double x, double y, double z, PreyType *pt) {
    /* This function gives an index of the number of pursuits per unit time made on items detected at each position.
     * However, it's based on probability densities, so it can't really translate into meaningful units unless it's
     * integrated over some volume. */
    if (z > surface_z || z < bottom_z) { return 0; }
    double prey_radius = pt->get_max_attended_distance();
    const double xsq = gsl_pow_2(x);
    const double ysq = gsl_pow_2(y);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(prey_radius);
    if (xsq + ysq + zsq > rsq) { return 0; }
    const double y0 = sqrt(rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    if (y < yT || y > y0) { return 0; }
    const double v = water_velocity(z);
    const double t_y = (y0 - y) / v; // Time at which the particle reaches position y
    auto integrand = [this, x, z, pt](double t)->double{
        double tauval = tau(t, x, z, pt);
        return exp(-t / tauval) / tauval;
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
            gsl_integration_qags(F, 0, t_y, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
            gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, 0, t_y, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    const double probability_of_no_previous_detection = 1 - result;
    const double tau_y = tau(t_y, x, z, pt);
    const double probability_of_detection_at_y = exp(-t_y / tau_y) / tau_y; // probability of detecting in one second (detection PDF) at position y
    const double answer = probability_of_no_previous_detection * probability_of_detection_at_y * v * pt->prey_drift_concentration * pt->true_hit_probability;
    return answer;
}

double Forager::relative_pursuits_by_position(double x, double y, double z) {
    /* Adds results from the above function across all prey types */
    double total = 0;
    for (auto &pt : prey_types) {
        total += relative_pursuits_by_position_single_prey_type(x, y, z, &pt);
    }
    return total;
}

std::map<std::string, std::vector<std::map<std::string, double>>> Forager::spatial_detection_proportions(PreyType *pt, std::string which_items, bool verbose) {
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
    std::vector<PreyType>* types;
    if (pt == nullptr) {
        types = &prey_types;                              // Run through all types if type was given as "All"
    } else {
        std::vector<PreyType> single_type_vector = {pt};   // Otherwise just do one, but make it a vector to use the same code
        types = &single_type_vector;
    }
    auto main_lambda = [this, which_items](PreyType *ipt, double min_distance, double max_distance, double min_angle, double max_angle) {  // integrates same function as above, but over the whole search volume
        double inner_theta, inner_phi;  // placeholders for the inner integrands
        double max_possible_distance = ipt->get_max_attended_distance();
        double truncated_min_distance = (min_distance > max_possible_distance) ? max_possible_distance : min_distance;
        double truncated_max_distance = (max_distance > max_possible_distance) ? max_possible_distance : max_distance;
        if (truncated_min_distance == truncated_max_distance) { return 0.0; }
        auto integrand = [this, ipt, which_items](double rho, double *theta, double *phi)->double{
            cartesian_3D_coords coords = cartesian_from_spherical(rho, *theta, *phi);
            if (coords.z > surface_z || coords.z < bottom_z) { return 0; }
            double v = water_velocity(coords.z);
            double prob = depletable_detection_probability(coords.x, coords.y, coords.z, ipt) * v;
            if (which_items == "Prey") {
                prob *= ipt->prey_drift_concentration * ipt->true_hit_probability;
            } else if (which_items == "Debris") {
                prob *= ipt->debris_drift_concentration * ipt->false_positive_probability;
            } else if (which_items == "All") {
                prob *= (ipt->prey_drift_concentration * ipt->true_hit_probability + ipt->debris_drift_concentration * ipt->false_positive_probability);
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
    for (auto &jpt : *types) {
        total += main_lambda(&jpt, NAN, NAN, NAN, NAN);
    }
//    std::vector<std::future<double>> total_futures;
//    for (auto &pt : *types) {
//        total_futures.push_back(std::async(std::launch::async, main_lambda, &pt, NAN, NAN, NAN, NAN));  // works fine, despite CLion highlighting error
//    };
//    for (auto &tf : total_futures) {
//        total += tf.get();
//    }
    std::vector<std::map<std::string, double>> distance_results, angle_results;
    for (auto distance_pair : distance_bins) {
        double min_distance = distance_pair.first;
        double max_distance = distance_pair.second;
        double bin_total = 0;
//        std::vector<std::future<double>> bin_futures;
//        for (auto &pt : *types) {
//            bin_futures.push_back(std::async(std::launch::async, main_lambda, &pt, min_distance, max_distance, NAN, NAN));  // works fine, despite CLion highlighting error
//        };
//        for (auto &bf : bin_futures) {
//            bin_total += bf.get();;
//        }
        for (auto &jpt : *types) {
            bin_total += main_lambda(&jpt, min_distance, max_distance, NAN, NAN);
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
//        std::vector<std::future<double>> bin_futures;
//        for (auto &pt : *types) {
//            bin_futures.push_back(std::async(std::launch::async, main_lambda, &pt, NAN, NAN, min_angle, max_angle));  // works fine, despite CLion highlighting error
//        };
//        for (auto &bf : bin_futures) {
//            bin_total += bf.get();;
//        }
        for (auto &jpt : *types) {
            bin_total += main_lambda(&jpt, NAN, NAN, min_angle, max_angle);
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

double Forager::proportion_of_detections_within(double min_distance, double max_distance, double min_angle, double max_angle,
                                                PreyType *pt, std::string *which_items) {
    /* Pass PreyType = NULL to get overall detections -- right now, that's the default. Want to allow specific categories.
     * Pass min_distance and max_distance = NAN to calculate proportions within a given angle.
     * Pass min_angle and max_angle = NAN to calculate proportions within a given distance. */

    //todo When I make this faster, add a verbose option that prints out the component proportions one-by-one. Calculate total first, so others are proportions.
    //todo Maybe parallelize, too? or will there be cache write conflicts / NaNs that require locks? shouldn't be many writes here...

    assert((isnan(min_distance) && isnan(max_distance)) || (min_distance >= 0 && max_distance <= max_radius));
    assert((isnan(min_angle) && isnan(max_angle)) || (min_angle >= 0 && max_angle <= theta));
    assert(*which_items == "Prey" || *which_items == "Debris" || *which_items == "All");
    if ((min_distance == max_distance) || (min_angle == max_angle)) {
        return 0;   // Sometimes comes up when trucating intended distance categories to fit within the foraging radius
    }
    auto region_lambda = [this, min_distance, max_distance, min_angle, max_angle, which_items](PreyType *ipt) { // integrates over the requested slice of the search volume
        double inner_theta, inner_phi;  // placeholders for the inner integrands
        auto integrand = [this, ipt, which_items](double rho, double *theta, double *phi)->double{
            cartesian_3D_coords coords = cartesian_from_spherical(rho, *theta, *phi);
            if (coords.z > surface_z || coords.z < bottom_z) { return 0; }
            double v = water_velocity(coords.z);
            double prob = depletable_detection_probability(coords.x, coords.y, coords.z, ipt) * v;
            if (*which_items == "Prey") {
                prob *= ipt->prey_drift_concentration * ipt->true_hit_probability;
            } else if (*which_items == "Debris") {
                prob *= ipt->debris_drift_concentration * ipt->false_positive_probability;
            } else if (*which_items == "All") {
                prob *= (ipt->prey_drift_concentration * ipt->true_hit_probability + ipt->debris_drift_concentration * ipt->false_positive_probability);
            }
            return prob;
        };
        gsl_function_pp_3d<decltype(integrand)> Fintegrand(integrand, &inner_theta, &inner_phi);
        auto *F = static_cast<gsl_function*>(&Fintegrand);
        double result = integrate_over_volume(F, min_distance, max_distance, min_angle, max_angle);
        if (isnan(result)) {
            printf("******** NaN returned for distance (%.8f,%.8f) for prey type %s items %s.\n", min_distance, max_distance, ipt->name.c_str(), which_items->c_str());
        }
        return result;
    };

    auto total_lambda = [this, which_items](PreyType *ipt) {  // integrates same function as above, but over the whole search volume
        double inner_theta, inner_phi;  // placeholders for the inner integrands
        auto integrand = [this, ipt, which_items](double rho, double *theta, double *phi)->double{
            cartesian_3D_coords coords = cartesian_from_spherical(rho, *theta, *phi);
            if (coords.z > surface_z || coords.z < bottom_z) { return 0; }
            double v = water_velocity(coords.z);
            double prob = depletable_detection_probability(coords.x, coords.y, coords.z, ipt) * v;
            if (*which_items == "Prey") {
                prob *= ipt->prey_drift_concentration * ipt->true_hit_probability;
            } else if (*which_items == "Debris") {
                prob *= ipt->debris_drift_concentration * ipt->false_positive_probability;
            } else if (*which_items == "All") {
                prob *= (ipt->prey_drift_concentration * ipt->true_hit_probability + ipt->debris_drift_concentration * ipt->false_positive_probability);
            }
            return prob;
        };
        gsl_function_pp_3d<decltype(integrand)> Fintegrand(integrand, &inner_theta, &inner_phi);
        auto *F = static_cast<gsl_function*>(&Fintegrand);
        return integrate_over_volume(F, NAN, NAN, NAN, NAN);
    };

    std::vector<PreyType> types;
    if (pt == nullptr) {
        types = prey_types;                              // Run through all types if type was given as "All"
    } else {
        std::vector<PreyType> single_type_vector = {pt};   // Otherwise just do one, but make it a vector to use the same code
        types = single_type_vector;
    }
    double region_result = 0;
    double total_result = 0;
    for (auto &ipt : types) {
        region_result += region_lambda(&ipt);
        total_result += total_lambda(&ipt);
    };
    return region_result / total_result;
}

PreyType *Forager::get_favorite_prey_type() {
    // Returns the favorite prey type by quantity of prey pursued
    PreyType *favorite;
    double favorite_pursuit_rate = -1;
    for (auto &pt : prey_types) {
        double rate = pursuit_rate("prey", &pt);
        if (rate > favorite_pursuit_rate) {
            favorite = &pt;
            favorite_pursuit_rate = rate;
        }
    }
    return favorite;
}

double Forager::pursuit_rate(std::string which_rate, PreyType *pt) {
    /* For which_rate, pass either "prey" or "debris." For prey type, pass a specific prey type object or nullptr
     * to sum over all prey types. */
    double x, y; // temporary holder for x values (y is required but unused) during the integrations
    assert(which_rate == "prey" || which_rate == "debris");
    std::vector<PreyType>* types;
    if (pt == nullptr) {
        types = &prey_types;                              // Run through all types if type was given as NULL
    } else {
        std::vector<PreyType> single_type_vector = {*pt};   // Otherwise just do one, but make it a vector to use the same code
        types = &single_type_vector;
    }
    auto numerator_integrand = [this, types, which_rate](double z, double *y, double *x) -> double {
        ++numerator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        double pd, pf, ph, dp, dd;
        for (auto &ipt : *types) {
            if (&ipt == nullptr) {
                printf("ERROR: Somehow inner prey type is null in numerator!!! Called with %lu types.\n", (*types).size());
            }
            pd = detection_probability(*x, z, &ipt);
            pf = ipt.false_positive_probability;
            ph = ipt.true_hit_probability;
            dp = ipt.prey_drift_concentration;
            dd = ipt.debris_drift_concentration;
            if (which_rate == "prey") {
                sum += (pd * v * (ph * dp));
            } else if (which_rate == "debris") {
                sum += (pd * v * (pf * dd));
            }
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(numerator_integrand)> Fn(numerator_integrand, &y, &x);
    gsl_function *Fnumerator = static_cast<gsl_function *>(&Fn);
    double numerator = integrate_over_xz_plane(Fnumerator, false);
    auto denominator_integrand = [this, types](double z, double *y, double *x) -> double {
        ++denominator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        double pd, pf, ph, h, dp, dd;
        for (auto &ipt : *types) {
            if (&ipt == nullptr) {
                printf("ERROR: Somehow inner prey type is null in denominator!!! Called with %lu types.\n", (*types).size());
            }
            pd = detection_probability(*x, z, &ipt);
            pf = ipt.false_positive_probability;
            ph = ipt.true_hit_probability;
            h = mean_maneuver_cost(*x, z, &ipt, false, pd);
            dp = ipt.prey_drift_concentration;
            dd = ipt.debris_drift_concentration;
            sum += pd * v * h * (ph * dp + pf * dd);
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
    if (DIAG_NANCHECKS) { assert(!isnan(rate)); }
    return rate;
}

void Forager::analyze_results() {
    /* Master command to compute all the metrics relevant for comparing the model to the data, following an optimization. */
    prey_pursuit_rate = pursuit_rate("prey", nullptr);
    debris_pursuit_rate = pursuit_rate("debris", nullptr);
    foraging_attempt_rate = prey_pursuit_rate + debris_pursuit_rate;
    proportion_of_attempts_ingested = prey_pursuit_rate / foraging_attempt_rate;
    for (auto &pt : prey_types) {
        pt.prey_pursuit_rate = pursuit_rate("prey", &pt);
        pt.debris_pursuit_rate = pursuit_rate("debris", &pt);
        pt.foraging_attempt_rate = pt.prey_pursuit_rate + pt.debris_pursuit_rate;
        pt.proportion_of_attempts_ingested = pt.prey_pursuit_rate / pt.foraging_attempt_rate;
        pt.diet_proportion = pt.prey_pursuit_rate / prey_pursuit_rate;
    }
}

