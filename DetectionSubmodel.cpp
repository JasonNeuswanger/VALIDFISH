//
// Created by Jason Neuswanger on 4/20/18.
//

#include "Forager.h"

inline double Forager::tau_effect_of_spatial_attention(double y, double distance) {
    const double halfFOV = theta / 2;
    const double angle = acos(y / distance);
    const double attention_dist = gsl_ran_gaussian_pdf(angle/sigma_A, sigma_A) / (sigma_A * (gsl_cdf_gaussian_P(halfFOV/sigma_A, sigma_A) - gsl_cdf_gaussian_P(-halfFOV/sigma_A, sigma_A)));
    const double effect = pow((1/theta)/attention_dist, A_0);
    assert(!isnan(effect));
    assert(effect >= 0);
    return effect;
}

inline double Forager::tau_effect_of_set_size() {
    // return pow(inspection_time * set_size, beta);
    const double effect = 1 + beta * set_size * inspection_time;
    assert(!isnan(effect));
    assert(effect >= 0);
}

inline double Forager::tau_effect_of_angular_area(double distance, const PreyType &pt) {
    // Calculate angular area, assuming the average prey item appears 4x as long as it is wide. The exact length:width
    // ratio doesn't really matter because it's scaled by a calibrated constant in the equation for tau anyway.
    const double angular_area = gsl_pow_2(atan2(pt.length, (M_PI * distance)));
    const double min_angular_area = 0.25 * gsl_pow_2(angular_resolution); // smallest visible
    if (angular_area <= min_angular_area) {
        printf("Angular area %.8f is <= min_angular_area %.8f, returning INF effect.\n", angular_area, min_angular_area);
        return INFINITY;
    } else {
        const double effect = delta_0 / (delta_0 + angular_area - min_angular_area);
        assert(!isnan(effect));
        assert(effect >= 0);
        return effect;
    }
}

inline double Forager::tau_effect_of_loom(double distance, double v, double y, const PreyType &pt) {
    // Calculate loom, the derivative of angular area with respect to time.
    const double loom = (2*M_PI*pt.length*v*y*atan(pt.length/(M_PI*distance))) / (distance * (gsl_pow_2(pt.length) + gsl_pow_2(M_PI * distance)));
    const double effect = (loom > 0) ? nu_0 / (nu_0 + loom) : 1;
    assert(!isnan(effect));
    assert(effect >= 0);
    return effect;
}

inline double Forager::tau_effect_of_search_image(const PreyType &pt) {
    const double effect = (pt.search_image_status == PreyType::SearchImageStatus::search_image_target) ? (1 / alpha_tau) : 1;
    assert(!isnan(effect));
    assert(effect >= 0);
    return effect;
}

double Forager::calculate_tau(double t, double x, double z, const PreyType &pt) {
    // Calculate some prerequisite quantities.
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = pt.rsq;
    const double v = water_velocity(z);
    assert(xsq + zsq <= rsq);
    if (xsq + zsq > rsq) { return INFINITY; }
    const double y = sqrt(rsq - xsq - zsq) - t * v;
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
    const double angular_length = 2 * atan2(pt.length, M_PI * distance);
    // We set tau=INFINITY to make the prey type & corresponding debris undetectable under the following circumstances:
    // 1) When angular size is too small to be visible at all.
    if (angular_length < angular_resolution) { return INFINITY; }
    // 2) When there is an active search image for a different prey type
    if (pt.search_image_status == PreyType::SearchImageStatus::search_image_exclusion) { return INFINITY; };
    // Combine the effects
    const double fish_attention_components = tau_effect_of_set_size() * tau_effect_of_spatial_attention(y, distance) * tau_effect_of_search_image(pt);
    const double object_salience_components = pt.crypticity * tau_effect_of_loom(distance, v, y, pt) * tau_effect_of_angular_area(distance, pt);
    const double combined_tau = (1/flicker_frequency) + tau_0 * fish_attention_components * object_salience_components;
    assert(!isnan(combined_tau));
    assert(combined_tau > 0);
    return combined_tau;
}

double Forager::tau(double t, double x, double z, const PreyType &pt) {
    double result;
    if (DIAG_NOCACHE or DIAG_NOCACHE_TAU) {
        result = calculate_tau(t, x, z, pt);
    } else {
        const double y = y_at_time(t, x, z, pt);
        long long key = xyzpt_hash_key(x, y, z, pt.uniqueid);
        auto cached_value = tau_cache.find(key);
        if (cached_value == tau_cache.end()) {
            result = calculate_tau(t, x, z, pt);
            tau_cache.insert(std::make_pair(key, result));
            ++tau_cache_misses;
        } else {
            ++tau_cache_hits;
            result = cached_value->second;
        }
    }
    return result;
}

std::map<std::string, double> Forager::tau_components(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    // Diagnostic version of the tau function to return the individual multipliers, used for plotting relative importance in Python.
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = pt->rsq;
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
    components.emplace("flicker_frequency", 1/flicker_frequency);
    components.emplace("tau_0", tau_0);
    components.emplace("crypticity", pt->crypticity);
    components.emplace("set_size", tau_effect_of_set_size());
    components.emplace("spatial_attention", tau_effect_of_spatial_attention(y, distance));
    components.emplace("search_image", tau_effect_of_search_image(*pt));
    components.emplace("loom", tau_effect_of_loom(distance, v, y, *pt));
    components.emplace("angular_area", tau_effect_of_angular_area(distance, *pt));
    return components;
};

void Forager::compute_set_size(bool verbose, char *printbuffer) {
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
            if (verbose && pt->prey_drift_concentration > 0.0) {
                sprintf(printbuffer+strlen(printbuffer), "For %20.20s, max. vis. dist=%.3f, set_volume=%.6f, prey_concentration=%4.1f, debris_concentration=%8.1f, ss for pt=%.3f.\n",
                       pt->name.c_str(), pt->max_visible_distance, set_volume, pt->prey_drift_concentration, pt->debris_drift_concentration, pt_ss);
            }
        }
    }
    set_size = ss;
}

double Forager::calculate_mean_value_function(double T, double x, double z, const PreyType &pt) {
    // This is the quantity referred to by a capital lambda (upside-down V) in the equations for a non-homogeneous Poisson
    // process. It is termed the "mean value function" in Sheldon's "Intro to Probability Models." It refers to the
    // expected value of the number of detections expected in the time interval (0, T). In this application, we
    // are only interested in the first detection, but this quantity is still used as an intermediate to calculate others.
    // Note that T here is not the full passage time, but the elapsed passage time to the position of interest, with T=0
    // being the time the object first becomes visible. It is integrated not necessarily from that point, but from the
    // point at which the item becomes profitable to pursue and therefore, by assumption, available for detection.
    // The upstream bound of integration is either the first time visible or the first time profitable, whichever is later,
    // and the lower bound is either the last time visible or the last time profitable, whichever is earlier.
    auto bounds = bounds_of_profitability(x, z, pt);
    if (isnan(bounds.first) || T <= bounds.first || bounds.first == bounds.second) { return 0; }
    const double upper_bound = (T > bounds.second) ? bounds.second : T;
    auto integrand = [this, x, z, pt](double t)->double{
        double tauval = tau(t, x, z, pt);
        return (isfinite(tauval)) ? 1/tauval : 0;
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
                gsl_integration_qags(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
                gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, bounds.first, upper_bound, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    assert(isfinite(result));
    assert(result >= 0);
    return result;
}

double Forager::mean_value_function(double T, double x, double z, const PreyType &pt) {
    double result;
    if (DIAG_NOCACHE || DIAG_NOCACHE_MEAN_VALUE_FUNCTION) {
        result = calculate_mean_value_function(T, x, z, pt);
    } else {
        const double y = y_at_time(T, x, z, pt);
        long long key = xyzpt_hash_key(x, y, z, pt.uniqueid);
        auto cached_value = mean_value_function_cache.find(key);
        if (cached_value == mean_value_function_cache.end()) {
            result = calculate_mean_value_function(T, x, z, pt);
            mean_value_function_cache.insert(std::make_pair(key, result));
            ++mean_value_function_cache_misses;
        } else {
            ++mean_value_function_cache_hits;
            result = cached_value->second;
        }
    }
    return result;
}

double Forager::calculate_detection_probability(double x, double z, const PreyType &pt) {
    const double T = passage_time(x, z, pt);
    const double result = 1 - exp(-mean_value_function(T, x, z, pt));
    assert(isfinite(result));
    assert(0 <= result && result <= 1);
    // ALTERNATIVE METHOD: Integrate the PDF. There is absolutely no reason to do this except for testing purposes, to
    // make sure the same result is obtained each way. Although that is the case, both results are subject to numerical
    // error reflected here. This can be arbitrarily reduced by decreasing QUAD_EPSREL at the expense of speed and/or
    // changing to adaptive quadrature. The result of that test suggests that we should use slightly higher precision
    // when integrating over the PDF to get the expected value of maneuver costs and discrimination probabilities than
    // during the other integrals.
//    auto integrand = [this, x, z, pt](double t)->double{
//        return detection_pdf_at_t(t, x, z, pt);
//    };
//    gsl_function_pp<decltype(integrand)> Fp(integrand);
//    gsl_function *F = static_cast<gsl_function*>(&Fp);
//    double result2, error;
//    #if USE_ADAPTIVE_INTEGRATION
//        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
//                    gsl_integration_qags(F, 0, T, QUAD_EPSABS, 0.1*QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result2, &error);
//                    gsl_integration_workspace_free(w);
//    #else
//        size_t neval;
//        gsl_integration_qng(F, 0, T, QUAD_EPSABS, 0.1*QUAD_EPSREL, &result2, &error, &neval);
//    #endif
//    assert(isfinite(result2));
//    assert(result2 >= 0);
//    if (result > 0) {
//        printf("Detection probability as direct CDF is %.8f. By integrating PDF, it's %.8f. Returning the former.",
//               result, result2);
//        if (abs(result - result2) > 1e-5) {
//            printf(" NOT A MATCH.");
//        }
//        printf("\n");
//    }
    return result;
}

double Forager::detection_probability(double x, double z, const PreyType &pt) {
    double result;
    if (DIAG_NOCACHE or DIAG_NOCACHE_DETECTION_PROBABILITY) {
        result = calculate_detection_probability(x, z, pt);
    } else {
        long long key = xzptiec_hash_key(x, z, pt.uniqueid, false);
        auto cached_value = detection_probability_cache.find(key);
        if (cached_value == detection_probability_cache.end()) {
            result = calculate_detection_probability(x, z, pt);
            detection_probability_cache.insert(std::make_pair(key, result));
            ++detection_probability_cache_misses;
        } else {
            ++detection_probability_cache_hits;
            result = cached_value->second;
        }
    }
    return result;
}

double Forager::calculate_detection_pdf_at_t(double t, double x, double z, const PreyType &pt) {
    if (z > surface_z || z < bottom_z) { return 0; }
    const double result = exp(-mean_value_function(t, x, z, pt)) / tau(t, x, z, pt);
    return result;
}

double Forager::detection_pdf_at_t(double t, double x, double z, const PreyType &pt) {
    double result;
    if (DIAG_NOCACHE or DIAG_NOCACHE_DETECTION_PDF) {
        result = calculate_detection_pdf_at_t(t, x, z, pt);
    } else {
        const double y = y_at_time(t, x, z, pt);
        long long key = xyzpt_hash_key(x, y, z, pt.uniqueid);
        auto cached_value = detection_pdf_cache.find(key);
        if (cached_value == detection_pdf_cache.end()) {
            result = calculate_detection_pdf_at_t(t, x, z, pt);
            detection_pdf_cache.insert(std::make_pair(key, result));
            ++detection_pdf_cache_misses;
        } else {
            ++detection_pdf_cache_hits;
            result = cached_value->second;
        }
    }
}

double Forager::detection_pdf_at_y(double y, double x, double z, const PreyType &pt) {
    const double t_y = time_at_y(y, x, z, pt);
    return detection_pdf_at_t(t_y, x, z, pt);
}

double Forager::detection_cdf_at_t(double t, double x, double z, const PreyType &pt) {    // Not used in model, just diagnostic
    return 1 - exp(-mean_value_function(t, x, z, pt));
}

double Forager::detection_cdf_at_y(double y, double x, double z, const PreyType &pt) {  // Not used in model, just diagnostic
    const double t_y = time_at_y(y, x, z, pt);
    return detection_cdf_at_t(t_y, x, z, pt);
}
