//
// Created by Jason Neuswanger on 4/20/18.
//

#include "Forager.h"

inline double Forager::tau_effect_of_spatial_attention(double y, double distance) {
    const double halfFOV = theta / 2;
    const double angle = acos(y / distance);
    const double attention_dist = gsl_ran_gaussian_pdf(angle/sigma_A, sigma_A) / (sigma_A * (gsl_cdf_gaussian_P(halfFOV/sigma_A, sigma_A) - (gsl_cdf_gaussian_P(-halfFOV/sigma_A, sigma_A))));
    double effect = pow((1/theta)/attention_dist, A_0);
    // printf("Tau effect of spatial attention is %.10f with y=%.3f, dist=%.3f, A_0=%.3f, angle=%.3f, term1=%.8f, term2=%.8f.\n", effect, y, distance, A_0, angle, gsl_cdf_gaussian_P(halfFOV/sigma_A, sigma_A), gsl_cdf_gaussian_P(-halfFOV/sigma_A, sigma_A));
    assert(isfinite(effect));
    assert(effect > 0);
    return effect;
}

inline double Forager::tau_effect_of_set_size() {
    return pow(inspection_time * set_size, beta);
}

inline double Forager::tau_effect_of_angular_area(double distance, std::shared_ptr<PreyType> pt) {
    // Calculate angular area, assuming the average prey item appears 4x as long as it is wide. The exact length:width
    // ratio doesn't really matter because it's scaled by a calibrated constant in the equation for tau anyway.
    const double angular_area = gsl_pow_2(atan2(pt->length, (M_PI * distance)));
    const double min_angular_area = 0.25 * gsl_pow_2(angular_resolution); // smallest visible
    if (angular_area <= min_angular_area) {
        return INFINITY;
    } else {
        return delta_0 / (delta_0 + angular_area - min_angular_area);
    }
}

inline double Forager::tau_effect_of_loom(double distance, double v, double y, std::shared_ptr<PreyType> pt) {
    // Calculate loom, the derivative of angular area with respect to time.
    const double loom = (2*M_PI*pt->length*v*y*atan(pt->length/(M_PI*distance))) / (distance * (gsl_pow_2(pt->length) + gsl_pow_2(M_PI * distance)));
    return (loom > 0) ? nu_0 / (nu_0 + loom) : 1;
}

inline double Forager::tau_effect_of_search_image(std::shared_ptr<PreyType> pt) {
    return (pt->search_image_status == PreyType::SearchImageStatus::search_image_target) ? (1 / alpha_tau) : 1;
}

inline double Forager::expected_profit_for_item(double t, double x, double y, double z, double v, std::shared_ptr<PreyType> pt){
    auto dps = discrimination_probabilities(t, x, z, pt);
    double false_positive_probability = dps.first;
    double true_hit_probability = dps.second;
    const double probability_a_pursued_item_is_prey = pt->prey_drift_concentration * true_hit_probability /
                                             (pt->prey_drift_concentration * true_hit_probability +
                                              pt->debris_drift_concentration * false_positive_probability);
    const double expected_energy_gain = pt->energy_content * probability_a_pursued_item_is_prey;
    const double maneuver_v = (v + focal_velocity) / 2;
    const double cost = maneuver_cost(x, y, z, maneuver_v, true);
    return expected_energy_gain - cost;
}

double Forager::tau(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    // Calculate some prerequisite quantities.
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(pt->get_max_visible_distance());
    const double v = water_velocity(z);
    if (xsq + zsq > rsq) {
        printf("TAU_ERROR: Somehow called for tau at (x,z)=(%.5f,%.5f) so xsq+rsq=%.5f is > rsq=%.5f. Prey type %s. Should be filtered out before here.\n", x, z, xsq+zsq, rsq, pt->name.c_str());
        return INFINITY;
    }
    const double y = sqrt(rsq - xsq - zsq) - t * v;
    if (isnan(y)) {
        printf("TAU_ERROR: y=nan at rsq=%.5f, xsq=%.5f, zsq=%.5f, t=%.5f, v=%.5f for prey type %s.\n", rsq, xsq, zsq, t, v, pt->name.c_str());
    }
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
    const double angular_length = 2 * atan2(pt->length, M_PI * distance);
    // We set tau=INFINITY to make the prey type & corresponding debris undetectable under the following circumstances:
    // 1) When angular size is too small to be visible at all.
    if (angular_length < angular_resolution) { return INFINITY; }
    // 2) When there is an active search image for a different prey type
    if (pt->search_image_status == PreyType::SearchImageStatus::search_image_exclusion) { return INFINITY; };
    // 3) When maneuvers would be expected to be unprofitable or impossible for the prey type at the given position. Obviously this
    //    does not actually render an item invisible (and therefore it is not excluded from set size calculations), but
    //    making it invisible here is the most computationally efficient/convenient place to exclude it from pursuit.
    if (exclude_unprofitable_maneuvers) {
        if (expected_profit_for_item(t, x, y, z, v, pt) <= 0) {
            return INFINITY;
        }
    } else {    // Even if we're not excluding unprofitable maneuvers, we still need to exclude impossible maneuvers, which return a cost of 100.
        const double maneuver_v = (v + focal_velocity) / 2;
        const double cost = maneuver_cost(x, y, z, maneuver_v, true);
        if (cost > 99) {
            return INFINITY;
        }
    }
    // Combine the effects
    const double fish_attention_components = tau_effect_of_set_size() * tau_effect_of_spatial_attention(y, distance) * tau_effect_of_search_image(pt);
    const double object_salience_components = pt->crypticity * tau_effect_of_loom(distance, v, y, pt) * tau_effect_of_angular_area(distance, pt);
    const double combined_tau = (1/flicker_frequency) + tau_0 * fish_attention_components * object_salience_components;
    assert(combined_tau > 0);
    if (combined_tau <= 0 or isnan(combined_tau)) {
        printf("TAU_ERROR: Combined_tau = %.3f. Attention part %.3f, salience %.3f. Components %.3f, %.3f, %.3f, %.3f, %.3f\n", combined_tau, fish_attention_components, object_salience_components, tau_effect_of_set_size(), tau_effect_of_spatial_attention(y, distance), tau_effect_of_search_image(pt), tau_effect_of_loom(distance, v, y, pt), tau_effect_of_angular_area(distance, pt));
    }
    return combined_tau;

}

std::map<std::string, double> Forager::tau_components(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    // Diagnostic version of the tau function to return the individual multipliers, used for plotting relative importance in Python.
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(pt->get_max_visible_distance());
    const double v = water_velocity(z);
    const double y = sqrt(rsq - xsq - zsq) - t * v;
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
    const double angular_length = 2 * atan2(pt->length, M_PI * distance);
    std::map<std::string, double> components;
    components.emplace("t", t); // not a component of tau, but an index for the others
    components.emplace("y", y); // not a component of tau, but gives context for the others
    const double standin_for_infinity = 10000;     // lower value to indicate infinite effect on the plots
    if (angular_length < angular_resolution) {
        components.emplace("angular_length_too_small_to_see", standin_for_infinity);
    } else {
        components.emplace("angular_length_too_small_to_see", 1);
    }
    if (pt->search_image_status == PreyType::SearchImageStatus::search_image_exclusion) {
        components.emplace("excluded_by_search_image", standin_for_infinity);
    } else {
        components.emplace("excluded_by_search_image", 1);
    }
    if (exclude_unprofitable_maneuvers && expected_profit_for_item(t, x, y, z, v, pt) <= 0) {
        components.emplace("excluded_for_unprofitabile_maneuver", standin_for_infinity);
    } else {
        components.emplace("excluded_for_unprofitable_maneuver", 1);
    }
    components.emplace("flicker_frequency", 1/flicker_frequency);
    components.emplace("tau_0", tau_0);
    components.emplace("crypticity", pt->crypticity);
    components.emplace("set_size", tau_effect_of_set_size());
    components.emplace("spatial_attention", tau_effect_of_spatial_attention(y, distance));
    components.emplace("search_image", tau_effect_of_search_image(pt));
    components.emplace("loom", tau_effect_of_loom(distance, v, y, pt));
    components.emplace("angular_area", tau_effect_of_angular_area(distance, pt));
    return components;
};

void Forager::compute_set_size(bool verbose) {
    // First, compute an adjustment for spatial attention that equals the integral under the portion of the attention
    // distribution PDF that lies below the level corresponding to a uniform attention distribution. In other words,
    // selectively unattended areas count less toward set size, but selectively more attended areas still only count
    // each item as 1 toward set size.
    const double halfFOV = theta / 2;
    const double attention_level_if_uniformly_distributed = 1 / theta;
    auto integrand = [this, halfFOV, attention_level_if_uniformly_distributed](double angle)->double{
        const double attention_dist = gsl_ran_gaussian_pdf(angle/sigma_A, sigma_A) / (sigma_A * (gsl_cdf_gaussian_P(halfFOV/sigma_A, sigma_A) - (gsl_cdf_gaussian_P(-halfFOV/sigma_A, sigma_A))));
        const double attention_effect_on_tau = pow(attention_level_if_uniformly_distributed/attention_dist, A_0);
        return fmin(1, 1/attention_effect_on_tau);
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    auto *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    size_t neval;
    gsl_integration_qng(F, -halfFOV, halfFOV, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    const double spatial_attention_adjustment = result / theta;
    // Now compute the rest of the set size
    double ss = 0;
    double set_volume, pt_ss;
    for (auto & pt : prey_types) {
        if (pt->search_image_status != PreyType::SearchImageStatus::search_image_exclusion) {
            set_volume = volume_within_radius(pt->max_visible_distance);
            pt_ss = set_volume * (pt->prey_drift_concentration + pt->debris_drift_concentration) * spatial_attention_adjustment;
            ss += pt_ss;
            if (verbose) {
                printf("For %20.20s, max. vis. dist=%.3f, set_volume=%.6f, prey_concentration=%4.1f, debris_concentration=%8.1f, ss for pt=%.3f.\n",
                       pt->name.c_str(), pt->max_visible_distance, set_volume, pt->prey_drift_concentration, pt->debris_drift_concentration, pt_ss);
            }
        }
    }
    set_size = ss;
}


double Forager::integrate_detection_pdf(double x, double z, std::shared_ptr<PreyType> pt) {
    /* Integrates the detection pdf over the prey's path from t=0 to t=T */
    printf("\n\nDoing an integration of the detection PDF at (x,z)=(%.3f,%.3f) for prey type %s:\n\n", x, z, pt->name.c_str());
    const double T = passage_time(x, z, pt);
    printf("Passage time is %.3f with water velocity %.5f.\n", T, water_velocity(z));
    if (T <= 0) { return 0; }  // If (x, z) are outside the search volume, detection probability is 0.
    auto integrand = [this, x, z, pt](double t)->double{
        double tauval = tau(t, x, z, pt);
        if (isfinite(tauval)) {
            printf("At t=%.3f, tau is %.5f. Detection PDF is %.5f.\n", t, tauval, exp(-t / tauval) / tauval);
            return exp(-t / tauval) / tauval;
        } else {
            printf("At t=%.3f, detection PDF is %.3f because tau is %.3f.\n", t, exp(-t / tauval) / tauval, tauval);
            return 0.0; // this is the value given by exp(-t / tauval) / tauval when tauval is inf
        }
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    // printf("Doing the integration.\n");
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
            gsl_integration_qags(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
            gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    assert(isfinite(result));
    printf("In integrate_detection_pdf, detection probability %.20f with error %.20f.\n", result, error);

//    if (!(0 <= result && result <= 1)) {  // todo figure out why some integrals here were > 1
//        printf("BAD result = %.5f.\n",result);
//        //abort();
//    } else {
//        printf("Good result = %.5f.\n",result);
//    }
    assert(0 <= result && result <= 1); // result here, when off, tends to range from 1.01 to 1.14 or so. Not improved by integration precision or adaptive algo though.
                                        // Problem only happens when excluding unprofitable (and/or impossible) maneuvers.
    return result;
}

double Forager::detection_probability(double x, double z, std::shared_ptr<PreyType> pt) {
    double result;
    if (DIAG_NOCACHE) {
        result = integrate_detection_pdf(x, z, pt);
    } else {
        long long key = xzpciec_hash_key(x, z, pt, false);
        auto cached_value = detection_probability_cache.find(key);
        if (cached_value == detection_probability_cache.end()) {
            result = integrate_detection_pdf(x, z, pt);
            detection_probability_cache[key] = result;
            ++detection_probability_cache_misses;
        } else {
            ++detection_probability_cache_hits;
            result = cached_value->second;
        }
    }
    if (DIAG_NANCHECKS) {
        assert(isfinite(result));
    }
    return result;
}