//
// Created by Jason Neuswanger on 5/12/18.
//

#include "Forager.h"

double Forager::integrate_energy_cost_over_prey_path(double x, double z, const PreyType &pt, bool is_energy_cost) {
    std::pair<double, double> bounds = bounds_of_profitability(x, z, pt);
    if (isnan(bounds.first) || bounds.first == bounds.second) { return 0; }
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    if (xsq + zsq >= pt.rhosq) { return 0; }    // if (x,z) are outside search volume return 0
    const double y0 = sqrt(pt.rsq - xsq - zsq);
    const double v = water_velocity(z); // Duplicating the passage_time function within this one because we need v below
    auto integrand = [this, x, z, pt, y0, v, is_energy_cost](double t)->double{
        const double y = y0 - v * t;
//        if (pt.rsq < gsl_pow_2(x) + gsl_pow_2(y) + gsl_pow_2(z)) {
//            printf("Bad one: %.5f + %.5f + %.5f = %.5f > %.5f.\n", gsl_pow_2(x), gsl_pow_2(y), gsl_pow_2(z), gsl_pow_2(x) + gsl_pow_2(y) + gsl_pow_2(z), pt.rsq);
//        }
//        assert(pt.rsq >= gsl_pow_2(x) + gsl_pow_2(y) + gsl_pow_2(z)); // Todo this assertion is sometimes violated by minor rounding errors, diagnose/fix
        const double det_pdf = detection_pdf_at_t(t, x, z, pt);
        if (det_pdf == 0) { return 0; }
        const double maneuver_v = (v + focal_velocity) / 2;  // Calculate cost from avg of focal & prey position velocities
        const double cost_if_pursued = maneuver_cost(x, y, z, maneuver_v, is_energy_cost);
        if (det_pdf * cost_if_pursued < 0) {
            printf("Got a detection PDF of %.12f with cost of %.12f leading to integrand of %.12f.\n", det_pdf, cost_if_pursued, det_pdf * cost_if_pursued);
        }
        return det_pdf * cost_if_pursued;
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
                gsl_integration_qags(F, 0, T, QUAD_EPSABS, 0.1*QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
                gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, bounds.first, bounds.second, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    assert(isfinite(result));
    if (result < 0) {
        printf("Got a negative energy cost, yikes: %.18f.\n", result);
    }
    assert(result >= 0);
    return result;
}

double Forager::expected_maneuver_cost(double x, double z, const PreyType &pt, bool is_energy_cost, double det_prob) {
    /* Right now this requires detection probability as an argument, because the only place it's called from is the
     * energy intake rate integral loop, and detection probability for the same x, z, pc is already precalculated
     * there. However, in the cached version at least, it's very fast to call the detection probability function
     * because it's already cached. If I need this function separately from that loop with precalculated detprob,
     * I can always put in a check to pass det_prob = -1 and force a recalculation if det_prob == -1. */
    double result;
    if (DIAG_NOCACHE) {
        result = (det_prob > 0.) ? integrate_energy_cost_over_prey_path(x, z, pt, is_energy_cost) / det_prob : 0.;
    } else {
        long long key = xzptiec_hash_key(x, z, pt.uniqueid, is_energy_cost);
        auto cached_value = expected_maneuver_cost_cache.find(key);
        if (cached_value == expected_maneuver_cost_cache.end()) {
            result = (det_prob > 0.) ? integrate_energy_cost_over_prey_path(x, z, pt, is_energy_cost) / det_prob : 0.;
            expected_maneuver_cost_cache.insert(std::make_pair(key, result));
            ++expected_maneuver_cost_cache_misses;
        } else {
            result = cached_value->second;
            ++expected_maneuver_cost_cache_hits;
        }
    }
    assert(isfinite(result));
    return result;
}

inline double Forager::item_profitability_at_time(double t, double x, double y, double z, const double maneuver_v, const PreyType &pt) {
    std::pair<double, double> dps = discrimination_probabilities(t, x, z, pt);
    const double false_positive_probability = dps.first;
    const double true_hit_probability = dps.second;
    const double probability_a_pursued_item_is_prey = pt.prey_drift_concentration * true_hit_probability /
                                                      (pt.prey_drift_concentration * true_hit_probability +
                                                       pt.debris_drift_concentration * false_positive_probability);
    const double expected_energy_gain = pt.energy_content * probability_a_pursued_item_is_prey;
    const double cost = maneuver_cost(x, y, z, maneuver_v, true);
    return expected_energy_gain - cost;
}

std::pair<double, double> Forager::calculate_bounds_of_profitability(double x, double z, const PreyType &pt) {
    /* Note that this method excludes the possibility that something would be profitable, then become unprofitable due
     * to something like higher angular velocity making identification harder, then become profitable again later on.
     * I think that's an unrealistic scenario for real fish, but I could see it potentially happening with the equations,
     * so I shuld try some sort of warning here if it occurs. */
    const double v = water_velocity(z);
    const double maneuver_v = (v + focal_velocity) / 2;
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    if (xsq + zsq > pt.rhosq) { std::make_pair(NAN, NAN); }    // If outside search volume, bounds are NAN
    const double y0 = sqrt(pt.rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    const double T = (y0 - yT) / v;
    double t, y, profit;
    double upstream_bound = NAN;
    double downstream_bound = NAN;
//    bool passed_downstream_bound = false;
    const int n_points = 50;
    bool downstream_bound_found = false;
    const double profit_at_downstream_visibility_limit = item_profitability_at_time(T, x, yT, z, maneuver_v, pt);
    if (profit_at_downstream_visibility_limit > 0) {    // We check if the downstream bound is the visibility limit
        downstream_bound = T;                           // To save the trouble of looping through all the other
        downstream_bound_found = true;                  // profitability calculations.
    }
    for (int i=0; i <= n_points; i++) {
        t = (i / (double) n_points) * T;
        y = y0 - t * v;
        profit = item_profitability_at_time(t, x, y, z, maneuver_v, pt);
        if (profit > 0) {
            if (isnan(upstream_bound)) {
                upstream_bound = t;
                if (downstream_bound_found) {
                    break;  // Exit the loop after finding the upstream bound if we already know the downstream bound is the visibility limit
                }
            }
            downstream_bound = t;
            // Todo This condition gets violated and prints these warnings when not commented. Figure out how extreme when using real parameters.
//            if (passed_downstream_bound) {
//                printf("Maneuver profitability bounds warning: maneuvers became profitable, then not, then profitable again farther downstream.\n");
//            }
        } else {
//            if (!isnan(downstream_bound)) {
//                passed_downstream_bound = true;
//            }
        }
    }
    return std::make_pair(upstream_bound, downstream_bound);
};

std::pair<double, double> Forager::bounds_of_profitability(double x, double z, const PreyType &pt) {
    if (DIAG_NOCACHE) {
        return calculate_bounds_of_profitability(x, z, pt);
    } else {
        long long key = xzptiec_hash_key(x, z, pt.uniqueid, false);
        auto cached_value = bounds_of_profitability_cache.find(key);
        if (cached_value == bounds_of_profitability_cache.end()) {
            std::pair<double, double> result = calculate_bounds_of_profitability(x, z, pt);
            bounds_of_profitability_cache.insert(std::make_pair(key, result));
            ++bounds_of_profitability_cache_misses;
            return result;
        } else {
            std::pair result = cached_value->second;
            ++bounds_of_profitability_cache_hits;
            return result;
        }
    }
}

bool Forager::location_is_profitable(double x, double y, double z, const PreyType &pt) {
    const double t = time_at_y(y, x, z, pt);
    auto bounds = bounds_of_profitability(x, z, pt);
    if (isfinite(bounds.first)) {
        return (t >= bounds.first && t <= bounds.second);
    } else {
        return false;
    }
}
