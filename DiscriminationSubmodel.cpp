//
// Created by Jason Neuswanger on 5/9/18.
//

#include "Forager.h"

inline double Forager::perception_effect_of_angular_area(double distance, const PreyType &pt) {
    const double angular_area = gsl_pow_2(atan2(pt.length, (M_PI * distance)));
    const double min_angular_area = 0.25 * gsl_pow_2(angular_resolution); // smallest visible
    if (angular_area <= min_angular_area) {
        return INFINITY;
    } else {
        return delta_p / (delta_p + angular_area - min_angular_area);
    }
}

inline double Forager::perception_effect_of_angular_velocity(double v, double t, double xsq, double zsq, double rsq) {
    const double angular_velocity = v * sqrt(xsq + xsq) / (rsq + t * v * (t * v - 2 * sqrt(rsq - xsq - zsq)));
    return (1 + pow(angular_velocity, omega_p)) / 2;
}

inline double Forager::perception_effect_of_search_image(const PreyType &pt) {
    double search_image_effect = 1;
    switch (pt.search_image_status) {
        case PreyType::SearchImageStatus::no_search_image:
            break; // stick with the default value of 1
        case PreyType::SearchImageStatus::search_image_exclusion:
            search_image_effect = 10000;
            break;
        case PreyType::SearchImageStatus::search_image_target:
            search_image_effect = alpha_d;
            break;
    }
    return search_image_effect;
}

inline double Forager::perception_effect_of_inspection_time() {
    return pow(1/inspection_time, ti_p);
}

double Forager::perceptual_variance(double t, double x, double z, const PreyType &pt) {
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double v = water_velocity(z);
    const double y = sqrt(pt.rsq - xsq - zsq) - t * v;
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
    return sigma_p_0
           * perception_effect_of_inspection_time()
           * perception_effect_of_search_image(pt)
           * perception_effect_of_angular_area(distance, pt)
           * perception_effect_of_angular_velocity(v, t, xsq, zsq, pt.rsq);
}

std::map<std::string, double> Forager::perceptual_variance_components(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    // Diagnostic version of the tau function to return the individual multipliers, used for plotting relative importance in Python.
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double v = water_velocity(z);
    const double y = sqrt(pt->rsq - xsq - zsq) - t * v;
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
    const double angular_length = 2 * atan2(pt->length, M_PI * distance);
    const double standin_for_infinity = 100;     // lower value to indicate infinite effect on the plots
    std::map<std::string, double> components;
    components.emplace("t", t); // not a component of perceptual variance, but an index for the others
    components.emplace("y", y); // not a component of perceptual variance, but gives context for the others
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
    components.emplace("sigma_p_0", sigma_p_0);
    components.emplace("inspection_time", perception_effect_of_inspection_time());
    components.emplace("angular_area", perception_effect_of_angular_area(distance, *pt));
    components.emplace("angular_velocity", perception_effect_of_angular_velocity(v, t, xsq, zsq, pt->rsq));
    return components;
};

std::pair<double, double> Forager::discrimination_probabilities(double t, double x, double z, const PreyType &pt) {
    // Saves some time over the function below by calculting both at once. But only worth doing when we can use both at once.
    if (DIAG_NO_DISCRIMINATION_MODEL) { return std::pair<double, double>{0.1, 0.9}; }
    if (DIAG_NOCACHE) {
        const double visual_sigma = sqrt(1 + perceptual_variance(t, x, z, pt));
        const double false_positive_probability = 1 - gsl_cdf_gaussian_P(discrimination_threshold / visual_sigma, 1);
        const double true_hit_probability =
                1 - gsl_cdf_gaussian_P((discrimination_threshold - discriminability) / visual_sigma, 1);
        return std::pair<double, double>{false_positive_probability, true_hit_probability};
    } else {
        const double y = y_at_time(t, x, z, pt);
        long long key = xyzpt_hash_key(x, y, z, pt.uniqueid);
        auto cached_value = discrimination_probability_cache.find(key);
        if (cached_value == discrimination_probability_cache.end()) {
            const double visual_sigma = sqrt(1 + perceptual_variance(t, x, z, pt));
            const double false_positive_probability =
                    1 - gsl_cdf_gaussian_P(discrimination_threshold / visual_sigma, 1);
            const double true_hit_probability =
                    1 - gsl_cdf_gaussian_P((discrimination_threshold - discriminability) / visual_sigma, 1);
            std::pair result = std::make_pair(false_positive_probability, true_hit_probability);
            discrimination_probability_cache.insert(std::make_pair(key, result));
            ++discrimination_probability_cache_misses;
            return result;
        } else {
            std::pair result = cached_value->second;
            ++discrimination_probability_cache_hits;
            return result;
        }
    }
}

double Forager::average_discrimination_probability_over_prey_path(double x, double z, const PreyType &pt, bool is_false_positive, double det_prob) {
    /* Finds the expected value of discrimination probability averaged over the prey's path from t=0 to t=T, given that
     * the item was eventaully detected. */
    auto bounds = bounds_of_profitability(x, z, pt);
    if (isnan(bounds.first) || bounds.first == bounds.second) { return 0; }
    auto integrand = [this, x, z, pt, is_false_positive](double t)->double{
        const double detection_pdf = detection_pdf_at_t(t, x, z, pt);
        const std::pair<double, double> disc_probs = discrimination_probabilities(t, x, z, pt);
        const double disc_prob = (is_false_positive) ? disc_probs.first : disc_probs.second;
        return disc_prob * detection_pdf;
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
    assert(0 <= result);
    // Capping the result at the detection probability corrects for cases in which minor (~1%) error in integrating the
    // pdf gives a value higher than the CDF, or the corresponding case for the expected value integral here.
    result = fmax(result, det_prob);
    return result / det_prob;
}

std::pair<double, double> Forager::expected_discrimination_probabilities(double x, double z, const PreyType &pt, double det_prob) {
    if (DIAG_NO_DISCRIMINATION_MODEL) { return std::pair<double, double> {0.1, 0.9}; }
    double false_positive_probability = 0;
    double true_hit_probability = 0;
    if (DIAG_NOCACHE) {
        if (det_prob > 0) {
            false_positive_probability = average_discrimination_probability_over_prey_path(x, z, pt, true, det_prob);
            true_hit_probability = average_discrimination_probability_over_prey_path(x, z, pt, false, det_prob);
            return std::pair<double, double>{false_positive_probability, true_hit_probability};
        } else {
            return std::pair<double, double>{0.5, 0.5};
        }
    } else {
        long long key = mdp_hash_key(x, z, pt.uniqueid, det_prob);
        auto cached_value = expected_discrimination_probability_cache.find(key);
        if (cached_value == expected_discrimination_probability_cache.end()) {
            if (det_prob > 0) {
                false_positive_probability = average_discrimination_probability_over_prey_path(x, z, pt, true, det_prob);
                true_hit_probability = average_discrimination_probability_over_prey_path(x, z, pt, false, det_prob);
            }
            std::pair result = std::make_pair(false_positive_probability, true_hit_probability);
            expected_discrimination_probability_cache.insert(std::make_pair(key, result));
            ++expected_discrimination_probability_cache_misses;
            return result;
        } else {
            std::pair result = cached_value->second;
            ++expected_discrimination_probability_cache_hits;
            return result;
        }
    }
}
