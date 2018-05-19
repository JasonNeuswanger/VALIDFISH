//
// Created by Jason Neuswanger on 1/31/18.
//

#include "Forager.h"

Forager::Forager(double fork_length_cm, double mass_g, double sigma_A, double mean_column_velocity,
                 double inspection_time, double discrimination_threshold, double search_image, double delta_0, double alpha_tau,
                 double alpha_d, double A_0, double beta, double bottom_z, double surface_z,
                 unsigned temperature_C, double tau_0, double flicker_frequency, double nu,
                 double discriminability, double delta_p, double omega_p, double ti_p, double sigma_p_0, std::string *maneuver_interpolation_csv_base_path)
        : Swimmer(fork_length_cm, mass_g, temperature_C, maneuver_interpolation_csv_base_path){

    #if GSL_ERROR_POLICY == 0
        gsl_set_error_handler(&ignore_gsl_errors); // Don't do anything about GSL convergence errors. Be careful with this in case of serious errors.
    #elif GSL_ERROR_POLICY == 1
        gsl_set_error_handler(&print_gsl_errors);  // Print an error instead of crashing when GSL had a convergence error.
    #else
        gsl_set_error_handler(NULL);               // Default to abort() on GSL errors. Useful for running the debugger.
    #endif
    compute_angular_resolution();
    set_strategy_bounds();
    set_parameter_bounds();
    // Fish strategy variables
    this->sigma_A = validate(s_sigma_A, sigma_A);
    this->mean_column_velocity = validate(s_mean_column_velocity, mean_column_velocity);
    this->inspection_time = validate(s_inspection_time, inspection_time);
    this->discrimination_threshold = validate(s_discrimination_threshold, discrimination_threshold);
    this->search_image = validate(s_search_image, search_image);
    // Model parameters that describe the fish's capabilities and need to be calibrated
    this->delta_0 = validate(p_delta_0, delta_0);
    this->alpha_d = validate(p_alpha_d, alpha_d);
    this->alpha_tau = validate(p_alpha_tau, alpha_tau);
    this->A_0 = validate(p_A_0, A_0);
    this->beta = validate(p_beta, beta);
    this->tau_0 = validate(p_tau_0, tau_0);
    this->flicker_frequency = validate(p_flicker_frequency, flicker_frequency);
    this->tau_0 = validate(p_tau_0, tau_0);
    this->nu_0 = validate(p_nu_0, nu);
    this->delta_p = validate(p_delta_p, delta_p);
    this->omega_p = validate(p_omega_p, omega_p);
    this->ti_p = validate(p_ti_p, ti_p);
    this->sigma_p_0 = validate(p_sigma_p_0, sigma_p_0);
    // Model parameters that describe the prey overall in unknown ways that need to be calibrated
    this->discriminability = validate(p_discriminability, discriminability);
    // Habitat variables
    this->bottom_z = bottom_z; assert(bottom_z < 0);             // river bottom z-coordinate in m (must be < 0)
    this->surface_z = surface_z; assert(surface_z > 0);          // river surface z-coordinate in m (must be > 0)
    this->depth = surface_z - bottom_z;
    // Physiological constraints from the literature, from from Wankowski (1979) as adapted by Hayes et al (2000)
    // Unlike other models, we don't truncated prey classes to suit these lengths; they're either excluded or not,
    // based on the mean size of prey in the category.
    this->min_prey_length_from_gill_rakers = 0.001 * (0.29 + 0.021 * this->fork_length_cm); // expressed here in meters
    this->max_prey_length_from_mouth_gape = 0.00105 * this->fork_length_cm * 4.3;    // same
}

Forager::Forager(Forager *otherForager) : Swimmer(*otherForager) {
    // Deep copy a forager and its prey categories and initialize new accelerators, caches, etc.
    sigma_A = otherForager->sigma_A;
    mean_column_velocity = otherForager->mean_column_velocity;
    inspection_time = otherForager->inspection_time;
    discrimination_threshold = otherForager->discrimination_threshold;
    // Model parameters
    tau_0 = otherForager->tau_0;
    delta_0 = otherForager->delta_0;
    alpha_tau = otherForager->alpha_tau;
    alpha_d = otherForager->alpha_d;
    A_0 = otherForager->A_0;
    beta = otherForager->beta;
    flicker_frequency = otherForager->flicker_frequency;
    nu_0 = otherForager->nu_0;
    discriminability = otherForager->discriminability;
    delta_p = otherForager->delta_p;
    omega_p = otherForager->omega_p;
    ti_p = otherForager->ti_p;
    sigma_p_0 = otherForager->sigma_p_0;
    // Habitat variables
    bottom_z = otherForager->bottom_z;
    surface_z = otherForager->surface_z;
    depth = surface_z - bottom_z;
    // Additional initialization
    for (auto & pt : otherForager->prey_types) {
        auto pt_copy = std::make_shared<PreyType>(&(*pt));
        prey_types.push_back(pt_copy); // creates a new prey category in-place in std::vector using copy constructor
    }
    angular_resolution = otherForager->angular_resolution;
    set_strategy_bounds();
    set_parameter_bounds();
    process_parameter_updates();
}

Forager::~Forager() {
    /* Destructor should deallocate any alloc'd c objects stored as instance variables, currently gsl interpolations */
}

void Forager::process_parameter_updates() {
    /* Carefully watch the ordering here, or else I could update some things with old values of other things. */
    Swimmer::process_parameter_updates(true);
    if (!DIAG_NOCACHE) {
        expected_maneuver_cost_cache.clear();
        detection_probability_cache.clear();
        expected_discrimination_probability_cache.clear();
        discrimination_probability_cache.clear();
        mean_value_function_cache.clear();
        tau_cache.clear();
        detection_pdf_cache.clear();
        bounds_of_profitability_cache.clear();
    }
    compute_angular_resolution();
    for (auto & pt : prey_types) {
        pt->compute_details(fork_length_cm, theta);
        double prey_type_radius = pt->max_visible_distance;
        if (prey_type_radius > max_radius) {
            max_radius = prey_type_radius;
        }
    }

    compute_focal_velocity();
    focal_swimming_cost = steady_swimming_cost(focal_velocity);
    compute_set_size(false);
}

void Forager::set_strategies(double sigma_A,
                             double mean_column_velocity,
                             double inspection_time,
                             double discrimination_threshold,
                             double search_image) {
    this->sigma_A = validate(s_sigma_A, sigma_A);
    this->mean_column_velocity = validate(s_mean_column_velocity, mean_column_velocity);
    this->inspection_time = validate(s_inspection_time, inspection_time);
    this->discrimination_threshold = validate(s_discrimination_threshold, discrimination_threshold);
    this->search_image = validate(s_search_image, search_image);
    process_parameter_updates();
}

void Forager::set_parameters(double delta_0,
                             double alpha_tau,
                             double alpha_d,
                             double beta,
                             double A_0,
                             double flicker_frequency,
                             double tau_0,
                             double nu_0,
                             double discriminability,
                             double delta_p,
                             double omega_p,
                             double ti_p,
                             double sigma_p_0) {
    this->delta_0 = delta_0;
    this->alpha_tau = alpha_tau;
    this->alpha_d = alpha_d;
    this->beta = beta;
    this->A_0 = A_0;
    this->tau_0 = tau_0;
    this->flicker_frequency = flicker_frequency;
    this->nu_0 = nu_0;
    this->discriminability = discriminability;
    this->delta_p = delta_p;
    this->omega_p = omega_p;
    this->ti_p = ti_p;
    this->sigma_p_0 = sigma_p_0;
    process_parameter_updates();
}

void Forager::set_strategy(Strategy strategy, double value) {
    switch (strategy) {
        case s_sigma_A: sigma_A = validate(s_sigma_A, value); break;
        case s_mean_column_velocity: mean_column_velocity = validate(s_mean_column_velocity, value); break;
        case s_inspection_time: inspection_time = validate(s_inspection_time, value); break;
        case s_discrimination_threshold: discrimination_threshold = validate(s_discrimination_threshold, value); break;
        case s_search_image: search_image = validate(s_search_image, value); break;
    }
    process_parameter_updates();
}

void Forager::set_parameter(Parameter parameter, double value) {
    switch (parameter) {
        case p_delta_0: delta_0 = validate(p_delta_0, value); break;
        case p_alpha_tau: alpha_tau = validate(p_alpha_tau, value); break;
        case p_alpha_d: alpha_d = validate(p_alpha_d, value); break;
        case p_beta: beta = validate(p_beta, value); break;
        case p_A_0: A_0 = validate(p_A_0, value); break;
        case p_flicker_frequency: flicker_frequency = validate(p_flicker_frequency, value); break;
        case p_tau_0: tau_0 = validate(p_tau_0, value); break;
        case p_nu_0: nu_0 = validate(p_nu_0, value); break;
        case p_discriminability: discriminability = validate(p_discriminability, value); break;
        case p_delta_p: delta_p = validate(p_delta_p, value); break;
        case p_omega_p: omega_p = validate(p_omega_p, value); break;
        case p_ti_p: ti_p = validate(p_ti_p, value); break;
        case p_sigma_p_0: sigma_p_0 = validate(p_sigma_p_0, value); break;
    }
    process_parameter_updates();
}

void Forager::compute_angular_resolution() {
    /* Computes the minimum angular size of prey a fish can detect, in radians. Based on using my angular size formula
     * and inverting the Hughes & Dill 1990 reaction distance equation (based on Schmidt & O'Brien 1992) to give angular
     * size. */
    angular_resolution = 2. * atan(1. / (120. * (1. - exp(-0.2 * fork_length_cm)) * M_PI));
}

double Forager::RateOfEnergyIntake(bool is_net, bool is_cost) {
    assert(num_prey_types() > 0);
    double x, y; // temporary holder for x values (y is required but unused) during the integrations
    auto numerator_integrand = [this, is_net, is_cost](double z, double *y, double *x)->double{
        ++numerator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        bool prey_is_resolvable;
        double pd, pf, ph, ch, E, dp, dd;
        for (auto & pt : prey_types) {
            prey_is_resolvable = (gsl_pow_2(*x) + gsl_pow_2(z) <= pt->rsq);
            if (prey_is_resolvable && (pt->prey_drift_concentration > 0 || pt->debris_drift_concentration > 0)) {
                pd = detection_probability(*x, z, *pt);
                auto mdps = expected_discrimination_probabilities(*x, z, *pt, pd);
                pf = mdps.first;     // false positive probability
                ph = mdps.second;    // true hit probability
                ch = (is_net || is_cost) ? expected_maneuver_cost(*x, z, *pt, true, pd) : 0;
                E = (is_cost) ? 0 : pt->energy_content;
                dp = pt->prey_drift_concentration;
                dd = pt->debris_drift_concentration;
                sum += (pd * v * ((E - ch) * ph * dp - ch * pf * dd));
            }
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(numerator_integrand)> Fn(numerator_integrand, &y, &x);
    auto Fnumerator = static_cast<gsl_function*>(&Fn);
    double net_energy_from_prey_maneuvers = integrate_over_xz_plane(Fnumerator, false);

    auto denominator_integrand = [this](double z, double *y, double *x)->double{
        ++denominator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        bool prey_is_resolvable;
        double pd, pf, ph, h, dp, dd;
        for (auto & pt : prey_types) {
            prey_is_resolvable = (gsl_pow_2(*x) + gsl_pow_2(z) <= pt->rsq);
            if (prey_is_resolvable && (pt->prey_drift_concentration > 0 || pt->debris_drift_concentration > 0)) {
                pd = detection_probability(*x, z, *pt);
                auto mdps = expected_discrimination_probabilities(*x, z, *pt, pd);
                pf = mdps.first;     // false positive probability
                ph = mdps.second;    // true hit probability
                // printf("For prey type %20s, pf=%.5f and ph=%.5f.\n", pt->name.c_str(), pf, ph);
                h = expected_maneuver_cost(*x, z, *pt, false, pd);
                dp = pt->prey_drift_concentration;
                dd = pt->debris_drift_concentration;
                sum += pd * v * h * (ph * dp + pf * dd);
            }
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(denominator_integrand)> Fd(denominator_integrand, &y, &x);
    auto Fdenominator = static_cast<gsl_function*>(&Fd);
    double total_handling_time = integrate_over_xz_plane(Fdenominator, false);
    double adjusted_focal_swimming_cost = (is_net) ? focal_swimming_cost : 0;
    double nrei = (net_energy_from_prey_maneuvers - adjusted_focal_swimming_cost) / (1 + total_handling_time);
    assert(!isnan(nrei));
    return nrei;
}

double Forager::NREI() {    // Net rate of energy intake
    return RateOfEnergyIntake(true, false);
}

double Forager::GREI() {    // Gross rate of energy intake
    return RateOfEnergyIntake(false, false);
}

double Forager::maneuver_cost_rate() {
    return -RateOfEnergyIntake(false, true); // negative sign changes the returned negative value to a positive cost value
}