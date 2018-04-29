//
// Created by Jason Neuswanger on 4/20/18.
//

#include "Forager.h"


// spatial attention has glitches.. getting y=nan, distance=nan, somehow!
// also had a glitch with angular are effect coming up negative

inline double Forager::tau_effect_of_spatial_attention(double y, double distance) {
    const double halfFOV = theta / 2;
    const double angle = acos(y / distance);
    const double attention_dist = gsl_ran_gaussian_pdf(angle/sigma_A, sigma_A) / (sigma_A * (gsl_cdf_gaussian_P(halfFOV/sigma_A, sigma_A) - (gsl_cdf_gaussian_P(-halfFOV/sigma_A, sigma_A))));
    const double effect = (A_0 + 1/theta) / (A_0 + attention_dist);
    if (isnan(effect)) {
        printf("TAU_ERROR: effect of spatial attentin is Nan with y=%.3f, dist=%.3f, A_0=%.3f, angle=%.3f, term1=%.8f, term2=%.8f.\n", y, distance, A_0, angle, gsl_cdf_gaussian_P(halfFOV/sigma_A, sigma_A), gsl_cdf_gaussian_P(-halfFOV/sigma_A, sigma_A));
    }
    return effect;
}

inline double Forager::tau_effect_of_set_size() {
    return 1 + beta * saccade_time * set_size;
}

inline double Forager::tau_effect_of_angular_area(double distance, std::shared_ptr<PreyType> pt) {
    // Calculate angular area, assuming the average prey item appears 4x as long as it is wide. The exact length:width
    // ratio doesn't really matter because it's scaled by a calibrated constant in the equation for tau anyway.
    const double angular_area = gsl_pow_2(atan2(pt->length, (M_PI * distance)));
    const double min_angular_area = 0.25 * gsl_pow_2(delta_min); // smallest visible
    if (angular_area <= min_angular_area) {
        return INFINITY;
    } else {
        return delta_0 / (angular_area - min_angular_area);
    }
}

inline double Forager::tau_effect_of_loom(double distance, double v, double y, std::shared_ptr<PreyType> pt) {
    // Calculate loom, the derivative of angular area with respect to time.
    const double loom = (2*M_PI*pt->length*v*y*atan(pt->length/(M_PI*distance))) / (distance * (gsl_pow_2(pt->length) + gsl_pow_2(M_PI * distance)));
    return (loom > 0) ? nu / (nu + loom) : 1;
}

inline double Forager::tau_effect_of_search_image(std::shared_ptr<PreyType> pt) {
    return (pt->search_image_status == PreyType::SearchImageStatus::search_image_target) ? (1 / alpha_tau) : 1;
}

double Forager::tau(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    // Calculate some prerequisite quantities.
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(pt->get_max_attended_distance());
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
    // 1) When angular size is too small to be visible at all, or too small for the current attention strategy (i.e.,
    //    angular_size < angular_resolution <= delta_min).
    if (angular_length < delta_min) { return INFINITY; }
    // 2) When there is an active search image for a different prey type
    if (pt->search_image_status == PreyType::SearchImageStatus::search_image_exclusion) { return INFINITY; };
    const double maneuver_v = (v + focal_velocity) / 2;
    // 3) When maneuvers would be expected to be unprofitable for the prey type at the given position. Obviously this
    //    does not actually render an item invisible (and therefore it is not excluded from set size calculations), but
    //    making it invisible here is the most computationally efficient/convenient place to exclude it from pursuit.
    if (exclude_unprofitable_maneuvers && maneuver_cost(x, y, z, maneuver_v, true) > pt->expected_energy_gain) { return INFINITY; };
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
    const double rsq = gsl_pow_2(pt->get_max_attended_distance());
    const double v = water_velocity(z);
    const double y = sqrt(rsq - xsq - zsq) - t * v;
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
    const double angular_length = 2 * atan2(pt->length, M_PI * distance);
    const double maneuver_v = (v + focal_velocity) / 2;
    std::map<std::string, double> components;
    components.emplace("t", t); // not a component of tau, but an index for the others
    components.emplace("y", y); // not a component of tau, but gives context for the others
    const double standin_for_infinity = 30;     // lower value to indicate infinite effect on the plots
    if (angular_length < delta_min) {
        components.emplace("angular_length_too_small_to_attend", standin_for_infinity);
    } else {
        components.emplace("angular_length_too_small_to_attend", 1);
    }
    if (pt->search_image_status == PreyType::SearchImageStatus::search_image_exclusion) {
        components.emplace("excluded_by_search_image", standin_for_infinity);
    } else {
        components.emplace("excluded_by_search_image", 1);
    }
    if (exclude_unprofitable_maneuvers && maneuver_cost(x, y, z, maneuver_v, true) > pt->expected_energy_gain) {
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