//
// Created by Jason Neuswanger on 2/8/18.
//

#include "Forager.h"

double Forager::depletable_detection_probability(double x, double y, double z, std::string *pc_name) {
    /* This function gives the probability of getting to a position undetected times the detection PDF at that
     * point. */
    if (z > surface_z || z < bottom_z) { return 0; }
    const double xsq = gsl_pow_2(x);
    const double ysq = gsl_pow_2(y);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(radius);
    if (xsq + ysq + zsq > rsq) { return 0; }
    const double y0 = sqrt(rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    if (y < yT || y > y0) { return 0; }
    const double v = water_velocity(z);
    const double t_y = (y0 - y) / v; // Time at which the particle reaches position y
    PreyCategory *pc = get_prey_category(pc_name);
    auto integrand = [this, x, z, pc](double t)->double{
        double tauval = tau(t, x, z, pc);
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
    const double tau_y = tau(t_y, x, z, pc);
    const double probability_of_detection_at_y = exp(-t_y / tau_y) / tau_y; // probability of detecting in one second (detection PDF) at position y
    const double answer = probability_of_no_previous_detection * probability_of_detection_at_y;
    // Commented line prints a ton of Python code to check that the Python model produces the same results
//    if (pc->uniqueid == 1) {
//        printf("print('Depletable detection probability for pc1 at (%.13f, %.13f, %.13f) is off by ', forager.depletable_detection_probability(%.13f, %.13f, %.13f, pc1) - %.13f)\n", x, y, z, x, y, z, answer);
//    }
    return answer;
}

double Forager::relative_pursuits_by_position_single_prey_type(double x, double y, double z, std::string *pc_name) {
    /* This function gives an index of the number of pursuits per unit time made on items detected at each position.
     * However, it's based on probability densities, so it can't really translate into meaningful units unless it's
     * integrated over some volume. */
    if (z > surface_z || z < bottom_z) { return 0; }
    const double xsq = gsl_pow_2(x);
    const double ysq = gsl_pow_2(y);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(radius);
    if (xsq + ysq + zsq > rsq) { return 0; }
    const double y0 = sqrt(rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    if (y < yT || y > y0) { return 0; }
    const double v = water_velocity(z);
    const double t_y = (y0 - y) / v; // Time at which the particle reaches position y
    PreyCategory *pc = get_prey_category(pc_name);
    auto integrand = [this, x, z, pc](double t)->double{
        double tauval = tau(t, x, z, pc);
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
    const double tau_y = tau(t_y, x, z, pc);
    const double probability_of_detection_at_y = exp(-t_y / tau_y) / tau_y; // probability of detecting in one second (detection PDF) at position y
    const double answer = probability_of_no_previous_detection * probability_of_detection_at_y * v * pc->prey_drift_density * pc->true_hit_probability;
    return answer;
}

double Forager::relative_pursuits_by_position(double x, double y, double z) {
    /* Adds results from the above function across all prey types */
    double total = 0;
    for (auto &pc : prey_categories) {
        total += relative_pursuits_by_position_single_prey_type(x, y, z, &pc.name);
    }
    return total;
}

double Forager::proportion_of_detections_within(double min_distance, double max_distance, double min_angle, double max_angle,
                                                std::string *pc_name, std::string *which_items) {
    /* Pass PreyCategory = NULL to get overall detections -- right now, that's the default. Want to allow specific categories.
     * Pass min_distance and max_distance = NAN to calculate proportions within a given angle.
     * Pass min_angle and max_angle = NAN to calculate proportions within a given distance. */
    assert((isnan(min_distance) && isnan(max_distance)) || (min_distance >= 0 && max_distance <= radius));
    assert((isnan(min_angle) && isnan(max_angle)) || (min_angle >= 0 && max_angle <= theta));
    assert(*which_items == "Prey" || *which_items == "Debris" || *which_items == "All");
    if ((min_distance == max_distance) || (min_angle == max_angle)) {
        return 0;   // Sometimes comes up when trucating intended distance categories to fit within the foraging radius
    }
    auto region_lambda = [this, min_distance, max_distance, min_angle, max_angle, which_items](PreyCategory *pc) { // integrates over the requeted slice of the search volume
        double inner_theta, inner_phi;  // placeholders for the inner integrands
        auto integrand = [this, pc, which_items](double rho, double *theta, double *phi)->double{
            cartesian_3D_coords coords = cartesian_from_spherical(rho, *theta, *phi);
            double v = water_velocity(coords.z);
            double prob = depletable_detection_probability(coords.x, coords.y, coords.z, &pc->name) * v;
            if (*which_items == "Prey") {
                prob *= pc->prey_drift_density * pc->true_hit_probability;
            } else if (*which_items == "Debris") {
                prob *= pc->debris_drift_density * pc->false_positive_probability;
            } else if (*which_items == "All") {
                prob *= (pc->prey_drift_density * pc->true_hit_probability + pc->debris_drift_density * pc->false_positive_probability);
            }
            return prob;
        };
        gsl_function_pp_3d<decltype(integrand)> Fintegrand(integrand, &inner_theta, &inner_phi);
        gsl_function *F = static_cast<gsl_function*>(&Fintegrand);
        return integrate_over_volume(F, min_distance, max_distance, min_angle, max_angle);
    };

    auto total_lambda = [this, which_items](PreyCategory *pc) {  // integrates same function as above, but over the whole search volume
        double inner_theta, inner_phi;  // placeholders for the inner integrands
        auto integrand = [this, pc, which_items](double rho, double *theta, double *phi)->double{
            cartesian_3D_coords coords = cartesian_from_spherical(rho, *theta, *phi);
            double v = water_velocity(coords.z);
            double prob = depletable_detection_probability(coords.x, coords.y, coords.z, &pc->name) * v;
            if (*which_items == "Prey") {
                prob *= pc->prey_drift_density * pc->true_hit_probability;
            } else if (*which_items == "Debris") {
                prob *= pc->debris_drift_density * pc->false_positive_probability;
            } else if (*which_items == "All") {
                prob *= (pc->prey_drift_density * pc->true_hit_probability + pc->debris_drift_density * pc->false_positive_probability);
            }
            return prob;
        };
        gsl_function_pp_3d<decltype(integrand)> Fintegrand(integrand, &inner_theta, &inner_phi);
        gsl_function *F = static_cast<gsl_function*>(&Fintegrand);
        return integrate_over_volume(F, NAN, NAN, NAN, NAN);
    };

    std::vector<PreyCategory>* categories;
    if (*pc_name == "All") {
        categories = &prey_categories;                              // Run through all categories if category was given as "All"
    } else {
        std::vector<PreyCategory> single_category_vector = {get_prey_category(pc_name)};   // Otherwise just do one, but make it a vector to use the same code
        categories = &single_category_vector;
    }
    double region_result = 0;
    double total_result = 0;
    for (auto &pc : *categories) {
        region_result += region_lambda(&pc);
        total_result += total_lambda(&pc);
    };
    // printf("Region total is %.16f and overall total is %.16f\n", region_result, total_result);
    return region_result / total_result;
}


double Forager::pursuit_rate(std::string which_rate, PreyCategory *pc) {
    /* For which_rate, pass either "prey" or "debris." For prey category, pass a specific category object or nullptr
     * to sum over all prey categories. */
    double x, y; // temporary holder for x values (y is required but unused) during the integrations
    assert(which_rate == "prey" || which_rate == "debris");
    std::vector<PreyCategory>* categories;
    if (pc == nullptr) {
        categories = &prey_categories;                              // Run through all categories if category was given as NULL
    } else {
        std::vector<PreyCategory> single_category_vector = {*pc};   // Otherwise just do one, but make it a vector to use the same code
        categories = &single_category_vector;
    }
    auto numerator_integrand = [this, categories, which_rate](double z, double *y, double *x) -> double {
        ++numerator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        double pd, pf, ph, dp, dd;
        for (auto &pc : *categories) {
            pd = detection_probability(*x, z, &pc);
            pf = pc.false_positive_probability;
            ph = pc.true_hit_probability;
            dp = pc.prey_drift_density;
            dd = pc.debris_drift_density;
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
    auto denominator_integrand = [this, categories](double z, double *y, double *x) -> double {
        ++denominator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        double pd, pf, ph, h, dp, dd;
        for (auto &pc : *categories) {
            pd = detection_probability(*x, z, &pc);
            pf = pc.false_positive_probability;
            ph = pc.true_hit_probability;
            h = mean_maneuver_cost(*x, z, &pc, false, pd);
            dp = pc.prey_drift_density;
            dd = pc.debris_drift_density;
            sum += pd * v * h * (ph * dp + pf * dd);
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(denominator_integrand)> Fd(denominator_integrand, &y, &x);
    gsl_function *Fdenominator = static_cast<gsl_function *>(&Fd);
    double denominator = integrate_over_xz_plane(Fdenominator, false);
    double rate = numerator / (1 + denominator);
    if (DIAG_NANCHECKS) { assert(!isnan(rate)); }
    return rate;
}

void Forager::analyze_results() {
    /* Master command to compute all the metrics relevant for comparing the model to the data, following an optimization. */
    prey_pursuit_rate = pursuit_rate("prey", nullptr);
    debris_pursuit_rate = pursuit_rate("debris", nullptr);
    foraging_attempt_rate = prey_pursuit_rate + debris_pursuit_rate;
    proportion_of_attempts_ingested = prey_pursuit_rate / foraging_attempt_rate;
    for (auto &pc : prey_categories) {
        pc.prey_pursuit_rate = pursuit_rate("prey", &pc);
        pc.debris_pursuit_rate = pursuit_rate("debris", &pc);
        pc.foraging_attempt_rate = pc.prey_pursuit_rate + pc.debris_pursuit_rate;
        pc.proportion_of_attempts_ingested = pc.prey_pursuit_rate / pc.foraging_attempt_rate;
        pc.diet_proportion = pc.prey_pursuit_rate / prey_pursuit_rate;
    }
}

