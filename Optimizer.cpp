//
// Created by Jason Neuswanger on 5/13/18.
//

#include "Optimizer.h"

Optimizer::Optimizer(std::shared_ptr<Forager> initial_forager, size_t max_iterations, size_t pack_size, bool verbose) {
    this->initial_forager = initial_forager;
    this->max_iterations = max_iterations;
    this->verbose = verbose;
    n_vars = (size_t) (1 + Forager::last_strategy - Forager::first_strategy);  // since the strategy enums are ints
    srand((unsigned) time(nullptr));  // Seed the random number generator for eigen todo might not be needed, check
    // Don't try to optimize search image if no prey types are eligible.
    if (initial_forager->num_search_image_eligible_prey_types() == 0) add_context(Forager::s_search_image, -1);
    // Create distributions of parameter bounds for random number generation.
    for (int sInt = Forager::first_strategy; sInt <= Forager::last_strategy; sInt++) {
        auto s = static_cast<Forager::Strategy>(sInt);
        distributions.insert({s, std::uniform_real_distribution<double> (initial_forager->strategy_bounds[s][0], initial_forager->strategy_bounds[s][1])});
    }
    // Create the wolves, which initialize themselves with random values.
    for (int i=0; i < pack_size; i++) {
        wolves.push_back(std::make_shared<Wolf>(this, &(*initial_forager), n_vars));
    }
    // Start the threads for all the wolves, which includes calculating their fitnesses.
    for (auto & wolf : wolves) {
        wolf->start_thread();
    }
    r1 = random_weights();
    r2 = random_weights();
}

double Optimizer::context_value(Forager::Strategy s) {
    if (context.find(s) != context.end()) {
        return context[s];
    } else {
        return NAN;
    }
}

double Optimizer::validated_random_strategy_value(Forager::Strategy s) {
    double context_val = context_value(s);
    if (isnan(context_val)) {
        return distributions[s](generator); // Generate a random number from a uniform distribution with the strategy's bounds
    } else {
        return context_val;
    }
}

void Optimizer::add_context(Forager::Strategy s, double value) {
    context[s] = value;
}

void Optimizer::clear_context() {
    // Clears all context settings and defaults to optimizing everything again, except it still won't optimize
    // search image if there aren't any prey types eligible for search images.
    context.clear();
    if (initial_forager->num_search_image_eligible_prey_types() == 0) add_context(Forager::s_search_image, -1);
}

ArrayXd Optimizer::random_weights() {
    return abs(ArrayXd::Random(n_vars, 1));
}

std::vector<double> Optimizer::optimize_forager() {
    std::vector<double> alpha_fitnesses_by_step;
    double a;     // Linearly decreases from 2 to 0 as max approaches
    int current_iteration = 0;
    wait_for_wolves_to_calculate_fitnesses();                   // Calculate fitnesses of initial wolves
    std::sort(wolves.begin(), wolves.end(),std::greater<>());
    alpha = &*wolves[0];
    beta = &*wolves[1];
    delta = &*wolves[2];
    printf("Beginning grey wolf optimization with initial alpha NREI of         %.8f J.\n", double(alpha->fitness));
    gsl_rng_default_seed = (unsigned long) time(nullptr);
    gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus);
    while (current_iteration < max_iterations) {
        a = 2. - 2. * ((double) current_iteration / max_iterations);
        // Create the weighting parameters for updating the next generation
        for (auto &wolf : wolves) {
            if (&*wolf != alpha && &*wolf != beta && &*wolf != delta) {
                r1 = random_weights();
                A1 = (2 * (a * r1)) - a;
                r2 = random_weights();
                C1 = 2 * r2;

                r1 = random_weights();
                A2 = (2 * (a * r1)) - a;
                r2 = random_weights();
                C2 = 2 * r2;

                r1 = random_weights();
                A3 = (2 * (a * r1)) - a;
                r2 = random_weights();
                C3 = 2 * r2;

                D_alpha = abs(C1 * alpha->params - wolf->params);    // Randomly weighted distance from the alpha
                D_beta = abs(C2 * beta->params - wolf->params);      // Randomly weighted distance from the beta
                D_delta = abs(C3 * delta->params - wolf->params);    // Randomly weighted distance from the delta
                X_1 = alpha->params - (A1 * D_alpha);
                X_2 = beta->params - (A2 * D_beta);
                X_3 = delta->params - (A3 * D_delta);
                wolf->set_params(0.5 * X_1 + 0.3 * X_2 + 0.2 * X_3);
            }
        }
        wait_for_wolves_to_calculate_fitnesses();
        // Note which wolves are the alpha, beta, and delta (top three) in this round.
        bool updated_alpha = false;
        if (current_iteration > 0) {
            for (auto &wolf : wolves) {
                if (*wolf > *alpha) {
                    delta = beta;
                    beta = alpha;
                    alpha = &*wolf;
                    updated_alpha = true;
                } else if (*wolf > *beta) {
                    delta = beta;
                    beta = &*wolf;
                } else if (*wolf > *delta) {
                    delta = &*wolf;
                };
            }
        }
        alpha_fitnesses_by_step.emplace_back(alpha->fitness);
        current_iteration += 1;
        if (verbose and updated_alpha) {
            printf("At the end of %4d iterations with %lu wolves, the alpha's fitness is %4.8f.          \n", current_iteration, wolves.size(), double(alpha->fitness));
        }
    }
    gsl_rng_free(rng);
    if (verbose) { printf("\n"); }
    // Update the initial forager object with the best solution
    initial_forager->set_strategies(
            alpha->params[Forager::Strategy::s_sigma_A],
            alpha->params[Forager::Strategy::s_mean_column_velocity],
            alpha->params[Forager::Strategy::s_inspection_time],
            alpha->params[Forager::Strategy::s_discrimination_threshold],
            alpha->params[Forager::Strategy::s_search_image]
    );
    // printf("Updating with alpha fitness %.5f with params %.5f, %.5f, %.5f, %.5f, %.5f.\n", double(alpha->fitness), alpha->params[0], alpha->params[1], alpha->params[2], alpha->params[3], alpha->params[4]);

    return alpha_fitnesses_by_step;
}

void Optimizer::wait_for_wolves_to_calculate_fitnesses() {
    bool wolves_are_done_computing_fitness = false;
    while (!wolves_are_done_computing_fitness) {
        wolves_are_done_computing_fitness = true;
        for (auto &wolf : wolves) {
            if (wolf->computing_fitness) {
                wolves_are_done_computing_fitness = false;
            }
        }
    }
}