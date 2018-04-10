//
// Created by Jason Neuswanger on 1/31/18.
//

#include "Forager.h"

Forager::Forager(double fork_length_cm,
                 double mass_g,
                 double radius,
                 double theta,
                 double mean_column_velocity,
                 double saccade_time,
                 double discrimination_threshold,
                 double delta_0,
                 double alpha_0,
                 double Z_0,
                 double c_1,
                 double beta,
                 double bottom_z,
                 double surface_z,
                 unsigned temperature_C,
                 double bed_roughness,
                 double discriminability,
                 double sigma_t,
                 double tau_0,
                 double t_V,
                 std::string *maneuver_interpolation_csv_base_path)
        : Swimmer(fork_length_cm, mass_g, temperature_C, maneuver_interpolation_csv_base_path){

    #if GSL_ERROR_POLICY == 0
        gsl_set_error_handler(&ignore_gsl_errors); // Don't do anything about GSL convergence errors. Be careful with this in case of serious errors.
    #elif GSL_ERROR_POLICY == 1
        gsl_set_error_handler(&print_gsl_errors);  // Print an error instead of crashing when GSL had a convergence error.
    #else
        gsl_set_error_handler(NULL);               // Default to abort() on GSL errors. Useful for running the debugger.
    #endif
    // Fish strategy variables
    this->radius = radius; assert(this->radius > 0);                  // radius of the the search volume
    this->theta = theta; assert(0 < this->theta < 2*M_PI);            // angular width of the search volume
    this->mean_column_velocity = mean_column_velocity; assert(this->mean_column_velocity > 0); // mean column velocity in m/s
    this->saccade_time = saccade_time;
    this->discrimination_threshold = discrimination_threshold;
    // Model parameters that describe the fish's capabilities and need to be calibrated
    this->delta_0 = delta_0;
    this->alpha_0 = alpha_0;
    this->Z_0 = Z_0;
    this->c_1 = c_1;
    this->beta = beta;
    this->t_V = t_V;
    // Model parameters that describe the prey overall in unknown ways that need to be calibrated
    this->discriminability = discriminability;
    this->sigma_t = sigma_t;
    this->tau_0 = tau_0;
    // Habitat variables
    this->bottom_z = bottom_z; assert(bottom_z < 0);             // river bottom z-coordinate in m (must be < 0)
    this->surface_z = surface_z; assert(surface_z > 0);          // river surface z-coordinate in m (must be > 0)
    this->bed_roughness = bed_roughness;                         // bed roughness height in m
    this->depth = surface_z - bottom_z;
    // Physiological constraints from the literature, from from Wankowski (1979) as adapted by Hayes et al (2000)
    // Unlike other models, we don't truncated prey classes to suit these lengths; they're either excluded or not,
    // based on the mean size of prey in the category.
    this->min_prey_length_from_gill_rakers = 0.000115 * this->fork_length_cm;        // converted to m
    this->max_prey_length_from_mouth_gape = 0.00105 * this->fork_length_cm * 4.3;    // same
    set_bounds();
}

Forager::Forager(Forager *otherForager) : Swimmer(*otherForager) {
    // Deep copy a forager and its prey categories and initialize new accelerators, caches, etc.
    radius = otherForager->radius;
    theta = otherForager->theta;
    mean_column_velocity = otherForager->mean_column_velocity;
    saccade_time = otherForager->saccade_time;
    discrimination_threshold = otherForager->discrimination_threshold;
    // Model parameters
    delta_0 = otherForager->delta_0;
    alpha_0 = otherForager->alpha_0;
    Z_0 = otherForager->Z_0;
    c_1 = otherForager->c_1;
    beta = otherForager->beta;
    t_V = otherForager->t_V;
    // Model parameters that describe the prey overall in unknown ways that need to be calibrated
    discriminability = otherForager->discriminability;
    sigma_t = otherForager->sigma_t;
    tau_0 = otherForager->tau_0;
    // Habitat variables
    bottom_z = otherForager->bottom_z;
    surface_z = otherForager->surface_z;
    bed_roughness = otherForager->bed_roughness;
    depth = surface_z - bottom_z;
    // Additional initialization
    for (auto &pc : otherForager->prey_categories) {
        //auto copied_prey_category = new PreyCategory(pc);
        //prey_categories.push_back(copied_prey_category);
        prey_categories.emplace_back(pc); // creates a new prey category in-place in std::vector using copy constructor
    }
    set_bounds();
    process_parameter_updates();
}

Forager::~Forager() {
    /* Destructor should deallocate any alloc'd c objects stored as instance variables, currently gsl interpolations */
};

void Forager::set_bounds() {
    bounds["radius"][0] = 0.01;
    bounds["radius"][1] = 10 * (0.01 * fork_length_cm); // ARBITRARY -- 10x fork length, converted to m
    bounds["theta"][0] = 0.01;
    bounds["theta"][1] = 1.85 * M_PI; // 1.85 pi is an estimate of fish's field of view; look for something in the literature
    bounds["mean_column_velocity"][0] = 0.01;
    bounds["mean_column_velocity"][1] = 10 * (0.01 * fork_length_cm);  // ARBITRARY
    bounds["saccade_time"][0] = 0.167; // physiological max 180 degrees/sec, assume 30 degree average, so 0.167 s just to move the eyes, let alone fixation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5050213/
    bounds["saccade_time"][1] = 1.0;   // just assuming the optimum will never be higher than this; adjust upward if needed
    bounds["discrimination_threshold"][0] = 0.;
    bounds["discrimination_threshold"][1] = 1.;
    bounds["attention"][0] = 0.;
    bounds["attention"][1] = 1.;
}

void Forager::modify_bound(std::string field, double lower_bound, double upper_bound) {
    bounds[field][0] = lower_bound;
    bounds[field][1] = upper_bound;
}

void Forager::fix_bound(std::string field, double fixed_value) {
    bounds[field][0] = fixed_value;
    bounds[field][1] = fixed_value;
}

void Forager::process_parameter_updates() {
    /* Carefully watch the ordering here, or else I could update some things with old values of other things. */
    Swimmer::process_parameter_updates(true);
    if (!DIAG_NOCACHE) {
        mean_maneuver_cost_cache.clear();
        detection_probability_cache.clear();
    }
    compute_angular_resolution();
    compute_discrimination_probabilities();
    compute_focal_velocity();
    focal_swimming_cost = steady_swimming_cost(focal_velocity);
    search_volume = volume_within_radius(radius);
    compute_search_rate();
    compute_set_size(false);
}

void Forager::modify_strategies(double radius,
                                double theta,
                                double mean_column_velocity,
                                double saccade_time,
                                double discrimination_threshold,
                                std::vector<double> attention) {

    assert(attention.size() == prey_categories.size());
    double total_attention = 0;
    for (auto &a : attention) { total_attention += a; }
    assert(fabs(total_attention - 1.) < 1e-12);   // Might not exactly equal 1 due to roundoff/floating point errors.
    this->radius = radius;
    this->theta = theta;
    this->mean_column_velocity = mean_column_velocity;
    this->saccade_time = saccade_time;
    this->discrimination_threshold = discrimination_threshold;
    for (size_t i=0; i < attention.size(); i++) {
        prey_categories.at(i).set_attention_allocated(attention.at(i));
    }
    process_parameter_updates();
}

void Forager::modify_parameters(double delta_0, double alpha_0, double beta, double Z_0, double c_1, double discriminability,
                                double sigma_t, double tau_0, double t_V) {
    this->delta_0 = delta_0;
    this->alpha_0 = alpha_0;
    this->beta = beta;
    this->Z_0 = Z_0;
    this->c_1 = c_1;
    this->discriminability = discriminability;
    this->sigma_t = sigma_t;
    this->tau_0 = tau_0;
    process_parameter_updates();
}

void Forager::modify_strategy(Strategy strategy, double value) {
    switch (strategy) {
        case s_radius: radius = value; break;
        case s_theta: theta = value; break;
        case s_mean_column_velocity: mean_column_velocity = value; break;
        case s_saccade_time: saccade_time = value; break;
        case s_discrimination_threshold: discrimination_threshold = value; break;
    }
    process_parameter_updates();
}

void Forager::modify_parameter(Parameter parameter, double value) {
    switch (parameter) {
        case p_delta_0: delta_0 = value; break;
        case p_alpha_0: alpha_0 = value; break;
        case p_beta: beta = value; break;
        case p_Z_0: Z_0 = value; break;
        case p_c_1: c_1 = value; break;
        case p_discriminability: discriminability = value; break;
        case p_sigma_t: sigma_t = value; break;
        case p_tau_0: tau_0 = value; break;
        case p_t_V: t_V = value; break;
    }
    process_parameter_updates();
}

double Forager::tau(double t, double x, double z, PreyCategory *pc) {
    double xsq = gsl_pow_2(x);
    double zsq = gsl_pow_2(z);
    double rsq = gsl_pow_2(radius);
    double v = water_velocity(z);
    double y = sqrt(rsq - xsq - zsq) - t * v;
    double angular_size = 2 * atan2(pc->length, (M_PI * sqrt(xsq + gsl_pow_2(y) + zsq)));
    if (angular_size < angular_resolution) { return INFINITY; }
    const double maneuver_v = (v + focal_velocity) / 2;
    if (exclude_unprofitable_maneuvers && maneuver_cost(x, y, z, maneuver_v, true) > pc->energy_content) { return INFINITY; }
    double retval = tau_0 * pc->crypticity * (1 + beta * saccade_time * set_size)
                          * (1 + alpha_0 / pc->get_feature_alpha())
                          * (1 + delta_0 / angular_size)
                          * (1 + search_rate / Z_0);
    assert(retval > 0);
    return retval;
}

double Forager::detection_probability(double x, double z, PreyCategory *pc) {
    double result;
    if (DIAG_NOCACHE) {
        result = integrate_detection_pdf(x, z, pc);
    } else {
        long long key = xzpciec_hash_key(x, z, pc, false);
        auto cached_value = detection_probability_cache.find(key);
        if (cached_value == detection_probability_cache.end()) {
            result = integrate_detection_pdf(x, z, pc);
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

void Forager::compute_angular_resolution() {
    /* Computes the minimum angular size of prey a fish can detect, in radians. Based on using my angular size formula
     * and inverting the Hughes & Dill 1990 reaction distance equation (based on Schmidt & O'Brien 1992) to give angular
     * size. */
    angular_resolution = 2. * atan(1. / (120. * (1. - exp(-0.2 * fork_length_cm)) * M_PI));
}

void Forager::compute_set_size(bool verbose) {
    double ss = 0;
    double set_alpha, set_radius, set_volume, pc_ss;
    for (auto & pc : prey_categories) {
        set_alpha = fmin(1, pc.get_feature_alpha()); // prevent any type from counting more than 1 per item toward set size
        if (set_alpha > 0) {
            set_radius = fmin(radius, pc.max_visible_distance(fork_length_cm)); // base set size on distance within which the item can be seen,
            set_volume = volume_within_radius(set_radius);                      // or the edge of the search radius, whichever is smaller
            pc_ss = set_alpha * set_volume * (pc.prey_drift_density + pc.debris_drift_density);
            ss += pc_ss;
            if (verbose) {
                printf("For %20.20s, max. vis. dist=%.3f, radius=%.3f, set_volume=%.6f, set_alpha=%.6f, prey_density=%4.1f, debris_density=%8.1f, ss for pc=%.3f.\n",
                    pc.name.c_str(), pc.max_visible_distance(fork_length_cm), set_radius, set_volume, set_alpha, pc.prey_drift_density, pc.debris_drift_density, pc_ss);
            }
        }
    }
    set_size = ss;
}

void Forager::compute_search_rate() {
    double volume_component = search_volume / t_V;
    auto velocity_xz = [this](double z)->double{ return water_velocity(z); };
    gsl_function_pp<decltype(velocity_xz)> Fp(velocity_xz);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double velocity_component = integrate_over_xz_plane(F, true);
    search_rate = volume_component + velocity_component;
}

void Forager::compute_discrimination_probabilities() {
    for (auto &pc : prey_categories) {
        double perceptual_sigma;
        if (pc.get_feature_alpha() <= 0 || saccade_time <= 0) {
            perceptual_sigma = 10000; // if alpha or ts are 0, make perceptual sigma huge (but not infinite, to avoid sqrt(0) below)
        } else {
            perceptual_sigma = sqrt(gsl_pow_2(sigma_t) + gsl_pow_2(c_1 / sqrt(pc.get_feature_alpha() * saccade_time)));
        }
        pc.false_positive_probability = 1 - gsl_cdf_gaussian_P(discrimination_threshold / perceptual_sigma, 1);
        pc.true_hit_probability = 1 - gsl_cdf_gaussian_P((discrimination_threshold - discriminability) / perceptual_sigma, 1);
    }
}

double Forager::mean_maneuver_cost(double x, double z, PreyCategory *pc, bool is_energy_cost, double det_prob) {
    /* Right now this requires detection probability as an argument, because the only place it's called from is the
     * energy intake rate integral loop, and detection probability for the same x, z, pc is already precalculated
     * there. However, in the cached version at least, it's very fast to call the detection probability function
     * because it's already cached. If I need this function separately from that loop with precalculated detprob,
     * I can always put in a check to pass det_prob = -1 and force a recalculation if det_prob == -1. */
    double result;
    if (DIAG_NOCACHE) {
        result = (det_prob > 0.) ? integrate_energy_cost_over_prey_path(x, z, pc, is_energy_cost) / det_prob : 0.;
    } else {
        long long key = xzpciec_hash_key(x, z, pc, is_energy_cost);
        auto cached_value = mean_maneuver_cost_cache.find(key);
        if (cached_value == mean_maneuver_cost_cache.end()) {
            result = (det_prob > 0.) ? integrate_energy_cost_over_prey_path(x, z, pc, is_energy_cost) / det_prob : 0.;
            mean_maneuver_cost_cache[key] = result;
            ++mean_maneuver_cost_cache_misses;
        } else {
            result = cached_value->second;
            ++mean_maneuver_cost_cache_hits;
        }
    }
    if (DIAG_NANCHECKS) { assert(isfinite(result)); }
    return result;
}

double Forager::RateOfEnergyIntake(bool is_net) {
    assert(num_prey_categories() > 0);
    double x, y; // temporary holder for x values (y is required but unused) during the integrations
    auto numerator_integrand = [this, is_net](double z, double *y, double *x)->double{
        ++numerator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        double pd, pf, ph, ch, E, dp, dd;
        for (auto & pc : prey_categories) {
            pd = detection_probability(*x, z, &pc);
            pf = pc.false_positive_probability;
            ph = pc.true_hit_probability;
            ch = (is_net) ? mean_maneuver_cost(*x, z, &pc, true, pd) : 0;
            E = pc.energy_content;
            dp = pc.prey_drift_density;
            dd = pc.debris_drift_density;
            sum += (pd * v * ((E - ch) * ph * dp - ch * pf * dd));
        }
        return sum;
    };
    gsl_function_pp_3d<decltype(numerator_integrand)> Fn(numerator_integrand, &y, &x);
    gsl_function *Fnumerator = static_cast<gsl_function*>(&Fn);
    double numerator = integrate_over_xz_plane(Fnumerator, false);

    auto denominator_integrand = [this](double z, double *y, double *x)->double{
        ++denominator_integrand_evaluations;
        double sum = 0;
        const double v = water_velocity(z);
        double pd, pf, ph, h, dp, dd;
        for (auto & pc : prey_categories) {
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
    gsl_function *Fdenominator = static_cast<gsl_function*>(&Fd);
    double denominator = integrate_over_xz_plane(Fdenominator, false);

    double swimming_cost = (is_net) ? focal_swimming_cost : 0;
    double nrei = (numerator - swimming_cost) / (1 + denominator);
    if (DIAG_NANCHECKS) { assert(!isnan(nrei)); }
    return nrei;
}

double Forager::NREI() {    // Net rate of energy intake
    return RateOfEnergyIntake(true);
}

double Forager::GREI() {    // Gross rate of energy intake
    return RateOfEnergyIntake(false);
}