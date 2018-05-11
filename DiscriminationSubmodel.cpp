//
// Created by Jason Neuswanger on 5/9/18.
//

#include "Forager.h"

inline double Forager::perception_effect_of_angular_area(double distance, std::shared_ptr<PreyType> pt) {
    const double angular_area = gsl_pow_2(atan2(pt->length, (M_PI * distance)));
    const double min_angular_area = 0.25 * gsl_pow_2(angular_resolution); // smallest visible
    if (angular_area <= min_angular_area) {
        return INFINITY;
    } else {
        return delta_p / (delta_p + angular_area - min_angular_area);
    }
}

inline double Forager::perception_effect_of_angular_velocity(double v, double t, double xsq, double zsq, double rsq, std::shared_ptr<PreyType> pt) {
    const double angular_velocity = v * sqrt(xsq + xsq) / (rsq + t * v * (t * v - 2 * sqrt(rsq - xsq - zsq)));
    return (1 + pow(angular_velocity, omega_p)) / 2;
}

inline double Forager::perception_effect_of_search_image(std::shared_ptr<PreyType> pt) {
    double search_image_effect = 1;
    switch (pt->search_image_status) {
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

double Forager::perceptual_variance(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(pt->get_max_visible_distance());
    const double v = water_velocity(z);
    const double y = sqrt(rsq - xsq - zsq) - t * v;
    const double ysq = gsl_pow_2(y);
    const double distance = sqrt(xsq + ysq + zsq);
//    printf("Perceptual variance = %.5f (from inspection time) * %.5f (from search image) * %.5f (from angular size) * %.5f (from angular velocity).\n", perception_effect_of_inspection_time(),
//           perception_effect_of_search_image(pt), perception_effect_of_angular_area(distance, pt), perception_effect_of_angular_velocity(v, t, xsq, zsq, rsq, pt));
    return sigma_p_0
           * perception_effect_of_inspection_time()
           * perception_effect_of_search_image(pt)
           * perception_effect_of_angular_area(distance, pt)
           * perception_effect_of_angular_velocity(v, t, xsq, zsq, rsq, pt);
}

std::pair<double, double> Forager::discrimination_probabilities(double t, double x, double z, std::shared_ptr<PreyType> pt) {
    // Saves some time over the function below by calculting both at once. But only worth doing when we can use both at once.
    const double visual_sigma = sqrt(1 + perceptual_variance(t, x, z, pt));
    const double false_positive_probability = 1 - gsl_cdf_gaussian_P(discrimination_threshold / visual_sigma, 1);
    const double true_hit_probability = 1 - gsl_cdf_gaussian_P((discrimination_threshold - discriminability) / visual_sigma, 1);
    // printf("For prey type %20s, at x=%.5f, z=%.5f, t=%.5f, pf=%.5f and ph=%.5f with prey concentration %.3f. Visual sigma=%.3f with perceptual variance = %.10f.\n", pt->name.c_str(), x, z, t, false_positive_probability, true_hit_probability, pt->prey_drift_concentration, visual_sigma, perceptual_variance(t, x, z, pt));
    return std::pair<double, double> { false_positive_probability, true_hit_probability };
}

double Forager::discrimination_probability(double t, double x, double z, std::shared_ptr<PreyType> pt, bool is_false_positive) {
    const double visual_sigma = sqrt(1 + perceptual_variance(t, x, z, pt));
    if (is_false_positive) {
        return 1 - gsl_cdf_gaussian_P(discrimination_threshold / visual_sigma, 1);
    } else  {
        return 1 - gsl_cdf_gaussian_P((discrimination_threshold - discriminability) / visual_sigma, 1);
    }
}

double Forager::average_discrimination_probability_over_prey_path(double x, double z, std::shared_ptr<PreyType> pt, bool is_false_positive, double det_prob) {
    /* Integrates the detection pdf over the prey's path from t=0 to t=T */
    const double T = passage_time(x, z, pt);
    if (T <= 0) { return 0.5; }  // If (x, z) are outside the search volume, discrimination is 50/50 luck.
    auto integrand = [this, x, z, pt, is_false_positive](double t)->double{
        const double tauval = tau(t, x, z, pt);
        const double detection_pdf = exp(-t / tauval) / tauval;
        const double disc_prob = discrimination_probability(t, x, z, pt, is_false_positive);
        return disc_prob * detection_pdf;
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
        gsl_integration_qng(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    printf("In average_discrimination_probability_over_prey_path, integration result %.5f, divide by detection probability %.5f. Used %lu evaluations.\n", result, det_prob, neval);
    assert(isfinite(result));
    //assert(0 <= result && result <= det_prob); // todo put this assert back in
    return result / det_prob;
}

std::pair<double, double> Forager::mean_discrimination_probabilities(double x, double z, std::shared_ptr<PreyType> pt, double det_prob) {
    double false_positive_probability = 0;
    double true_hit_probability = 0;
    if (true) {
        if (det_prob > 0) {
            assert(det_prob <= 1);
            false_positive_probability = average_discrimination_probability_over_prey_path(x, z, pt, true, det_prob);
            true_hit_probability = average_discrimination_probability_over_prey_path(x, z, pt, false, det_prob);
            //printf("For prey type %20s, pf=%.5f and ph=%.5f with prey drift concentration %.5f.\n", pt->name.c_str(), false_positive_probability, true_hit_probability, pt->prey_drift_concentration);
            return std::pair<double, double>{false_positive_probability, true_hit_probability};
        } else {
            return std::pair<double, double>{0.5, 0.5};
        }
    } else { // todo build cache for this after testing that things basically work
//        long long key = mdp_hash_key(x, z, pt, det_prob);
//        auto cached_value = mean_discrimination_probability_cache.find(key);
//        if (cached_value == mean_discrimination_probability_cache.end()) {
//            if (det_prob > 0) {
//                false_positive_probability = average_discrimination_probability_over_prey_path(x, z, pt, true, det_prob) / det_prob;
//                true_hit_probability = average_discrimination_probability_over_prey_path(x, z, pt, false, det_prob) / det_prob;
//            }
//            mean_discrimination_probability_cache[key] = std::pair<double, double> {false_positive_probability, true_hit_probability};
//            ++mean_discrimination_probability_cache_misses;
//        } else {
//            result = cached_value->second;
//            ++mean_discrimination_probability_cache_hits;
//        }
    }
}
