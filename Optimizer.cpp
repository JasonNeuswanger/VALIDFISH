//
// Created by Jason Neuswanger on 2/12/18.
//

#include "Optimizer.h"

Optimizer::Optimizer(Forager *initial_forager, size_t max_iterations, size_t pack_size, bool verbose) {
    this->initial_forager = initial_forager;
    this->max_iterations = max_iterations;
    this->pack_size = pack_size; // Probably ideally a multiple of the # of processor cores?
    this->verbose = verbose;
    n_vars = 6;
    srand(time(NULL));  // Seed the random number generator for the creation of random wolves
    if (initial_forager->num_search_image_eligible_prey_types() == 0) {
        add_context(Forager::s_search_image, -1);   // don't try to optimize search image if no prey types are eligible
    }
    for (int i=0; i < pack_size; i++) {
        wolves.push_back(random_wolf());    // Create the initial pack of wolves
    }
    r1 = random_weights();
    r2 = random_weights();
    calculate_wolf_fitnesses();
}

void Optimizer::set_algorithm_options(bool use_chaos, bool use_dynamic_C, bool use_exponential_decay, bool use_levy,
                                      bool use_only_alpha, bool use_weighted_alpha) {
    algorithm_use_chaos = use_chaos;                            // Chaotic Gray Wolf Optimization (Kohli & Arora 2017)
    algorithm_use_dynamic_C = use_dynamic_C;                    // Improved Grey Wolf Algorithm (Kumar et al 2016)
    algorithm_use_exponential_decay = use_exponential_decay;    // Modified Grey Wolf Optimizer (Mittal et al 2016)
    algorithm_use_levy = use_levy;                              // Levy Flight enhancement for GWO (Luo et al 2017)
    algorithm_use_only_alpha = use_only_alpha;                  // Enhanced Grey Wolf Optimization (Joshi & Arora 2017)
    algorithm_use_weighted_alpha = use_weighted_alpha;          // Custom variation of the above
}

double Optimizer::validated_random_parameter_value(Forager::Strategy s) {
    if (context.find(s) != context.end()) {
        return context[s];
    } else {
        return fRand(initial_forager->strategy_bounds[s][0], initial_forager->strategy_bounds[s][1]);
    }
}

void Optimizer::add_context(Forager::Strategy s, double value) {
    context[s] = value;
}

void Optimizer::clear_context() {
    // Clears all context settings and defaults to optimizing everything again, except it still won't optimize
    // search image if there aren't any prey types eligible for search images.
    context.clear();
    if (initial_forager->num_search_image_eligible_prey_types() == 0) {
        add_context(Forager::s_search_image, -1);   // don't try to optimize search image if no prey types are eligible
    }
}

wolf_type Optimizer::random_wolf() {    // Only used on initialization
    wolf_type wolf;
    wolf.params = ArrayXd::Random(n_vars, 1).abs(); // is this line necessary?
    wolf.params[0] = validated_random_parameter_value(Forager::s_delta_min);
    wolf.params[1] = validated_random_parameter_value(Forager::s_sigma_A);
    wolf.params[2] = validated_random_parameter_value(Forager::s_mean_column_velocity);
    wolf.params[3] = validated_random_parameter_value(Forager::s_saccade_time);
    wolf.params[4] = validated_random_parameter_value(Forager::s_discrimination_threshold);
    wolf.params[5] = validated_random_parameter_value(Forager::s_search_image);
    return wolf;
}

void Optimizer::enforce_bounds_and_constraints() {
    for (auto &wolf : wolves) {
        wolf.params[0] = trim_to_bounds(wolf.params[0], initial_forager->strategy_bounds[Forager::s_delta_min]);
        wolf.params[1] = trim_to_bounds(wolf.params[1], initial_forager->strategy_bounds[Forager::s_sigma_A]);
        wolf.params[2] = trim_to_bounds(wolf.params[2], initial_forager->strategy_bounds[Forager::s_mean_column_velocity]);
        wolf.params[3] = trim_to_bounds(wolf.params[3], initial_forager->strategy_bounds[Forager::s_saccade_time]);
        wolf.params[4] = trim_to_bounds(wolf.params[4], initial_forager->strategy_bounds[Forager::s_discrimination_threshold]);
        wolf.params[5] = trim_to_bounds(wolf.params[5], initial_forager->strategy_bounds[Forager::s_search_image]);
    }
}

ArrayXd Optimizer::random_weights() {
    return abs(ArrayXd::Random(n_vars, 1));
}

ArrayXd Optimizer::new_weights(ArrayXd previous_weights) {
    /* Keeps track of two random weights via chaotic maps, r1 and r2 */
    if (algorithm_use_chaos) {
        ArrayXd chaotic_arr(n_vars);
        for (int i = 0; i < n_vars; i++) {
            chaotic_arr[i] = (1 + cos(2 * acos(2 * previous_weights[i] - 1))) / 2.;
            assert(!isnan(chaotic_arr[i]));
            assert(chaotic_arr[i] >= 0);
            assert(chaotic_arr[i] <= 1);
        }
        return chaotic_arr;
    } else {
        return random_weights();
    }
}

std::vector<double> Optimizer::optimize_forager() {
    std::vector<double> alpha_fitnesses_by_step;
    double a;     // Linearly decreases from 2 to 0 as max approaches
    double a2;    // Modifications for another algorithm
    int current_iteration = 0;
    std::sort(wolves.begin(), wolves.end(),std::greater<>());
    alpha = &wolves[0];
    beta = &wolves[1];
    delta = &wolves[2];
    printf("Beginning grey wolf optimization with initial alpha NREI of         %.8f J.\n", alpha->fitness);
    const double exponent = (algorithm_use_exponential_decay) ? 2. : 1.;
    gsl_rng_default_seed = (unsigned long) time(nullptr);
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);

    while (current_iteration < max_iterations) {
        a = 2. - 2. * pow(((double) current_iteration / max_iterations), exponent);
        a2 = (algorithm_use_dynamic_C) ? a / 2. : 0.;
        // Create the weighting parameters for updating the next generation

        for (auto &wolf : wolves) {
            if (&wolf != alpha && &wolf != beta && &wolf != delta) {

                r1 = new_weights(r1);
                A1 = (2 * (a * r1)) - a;
                r2 = new_weights(r2);
                C1 = 2 * r2 - a2;

                r1 = new_weights(r1);
                A2 = (2 * (a * r1)) - a;
                r2 = new_weights(r2);
                C2 = 2 * r2 - a2;

                r1 = new_weights(r1);
                A3 = (2 * (a * r1)) - a;
                r2 = new_weights(r2);
                C3 = 2 * r2 - a2;

                D_alpha = abs(C1 * alpha->params - wolf.params);    // Randomly weighted distance from the alpha
                D_beta = abs(C2 * beta->params - wolf.params);      // Randomly weighted distance from the beta
                D_delta = abs(C3 * delta->params - wolf.params);    // Randomly weighted distance from the delta
                X_1 = alpha->params - (A1 * D_alpha);
                X_2 = beta->params - (A2 * D_beta);
                X_3 = delta->params - (A3 * D_delta);
                if (algorithm_use_only_alpha) {
                    wolf.params = X_1;
                } else if (algorithm_use_weighted_alpha) {
                    wolf.params = 0.5 * X_1 + 0.3 * X_2 + 0.2 * X_3;
                } else {
                    wolf.params = (X_1 + X_2 + X_3) / 3.;   // standard GWO equally weights all 3 top wolves
                }
                if (algorithm_use_levy) {    // previous levy method tests use scale 0.2, alpha 1.8
                    double levy_scale = 0.5;
                    double levy_alpha = 0.8; // can range from 0 < alpha <= 2; lower values = wider tails = more variation
                    for (int i=0; i < n_vars; i++) {
                        wolf.params[i] = wolf.params[i] + wolf.params[i] * gsl_ran_levy(rng, levy_scale, levy_alpha);
                    }

                }

            }
        }
        enforce_bounds_and_constraints();
        calculate_wolf_fitnesses();

        // Note which wolves are the alpha, beta, and delta (top three) in this round.
        bool updated_alpha = false;
        if (current_iteration > 0) {
            for (auto &wolf : wolves) {
                if (wolf.fitness > alpha->fitness) {
                    delta = beta;
                    beta = alpha;
                    alpha = &wolf;
                    updated_alpha = true;
                } else if (wolf.fitness > beta->fitness) {
                    delta = beta;
                    beta = &wolf;
                } else if (wolf.fitness > delta->fitness) {
                    delta = &wolf;
                };
            }
        }
        alpha_fitnesses_by_step.emplace_back(alpha->fitness);
        current_iteration += 1;
        if (verbose and updated_alpha) {
            printf("At the end of %4d iterations with %lu wolves, the alpha's fitness is %4.8f.          \n", current_iteration, wolves.size(), alpha->fitness);
        }
    }
    gsl_rng_free(rng);
    if (verbose) { printf("\n"); }
    // Update the initial forager object with the best solution
    update_forager_from_wolf(initial_forager, alpha);
    return alpha_fitnesses_by_step;
}

void Optimizer::calculate_wolf_fitnesses() {
    struct fitness_results { wolf_type * wolf; double fitness; };
    std::vector<std::future<fitness_results>> fitness_futures;
    auto wolf_fitness_lambda = [this](wolf_type *wolf) {
        struct fitness_results fr = {wolf, wolf_fitness(wolf)};
        return fr;
    };
    for (auto &wolf : wolves) {
        fitness_futures.push_back(std::async(std::launch::async, wolf_fitness_lambda, &wolf));
    };
    for (auto &result : fitness_futures) {
        fitness_results r = result.get();
        r.wolf->fitness = r.fitness;    // CLion IDE says this value is never used; yes it is!
    }
}

void Optimizer::update_forager_from_wolf(Forager *forager, wolf_type *wolf) {
    forager->set_strategies(
            wolf->params[0],       // delta_min
            wolf->params[1],       // sigma_A
            wolf->params[2],       // mean column velocity
            wolf->params[3],       // saccade time
            wolf->params[4],       // discrimination threshold
            wolf->params[5]        // search image
    );
}

double Optimizer::wolf_fitness(wolf_type *wolf) {
    /* Converts the wolf (an array of parameters) into a Forager object with the same parameters and calculates its
     * fitness (NREI). */
    auto forager = new Forager(initial_forager);  // Deep-copy a forager from the initial forager -- CONSIDER JUST UPDATING
    update_forager_from_wolf(forager, wolf);
    const double nrei = forager->NREI();
    delete forager;
    return nrei;
}