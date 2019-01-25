//
// Created by Jason Neuswanger on 5/14/18.
//

#include "Wolf.h"

Wolf::Wolf(Optimizer *optimizer, Forager *forager, size_t n_vars) : optimizer(optimizer), forager(Forager(forager)), the_thread() {
    params = Eigen::ArrayXd(n_vars);
    for (int sInt = Forager::first_strategy; sInt <= Forager::last_strategy; sInt++) {
        auto strategy = static_cast<Forager::Strategy>(sInt);
        params[strategy] = optimizer->validated_random_strategy_value(strategy);
    }
}

Wolf::~Wolf() {
    thread_running = false;
    if(the_thread.joinable()) the_thread.join();
}

void Wolf::enforce_bounds_and_context() {
    for (int sInt = Forager::first_strategy; sInt <= Forager::last_strategy; sInt++) {
        auto strategy = static_cast<Forager::Strategy>(sInt);
        const double context_val = optimizer->context_value(strategy);
        if (isnan(context_val)) {
            params[strategy] = trim_to_bounds(params[strategy], optimizer->initial_forager->strategy_bounds[strategy]);
        } else {
            params[strategy] = context_val;
        }
    }
}

void Wolf::calculate_fitness() {
    computing_fitness = true;
    enforce_bounds_and_context();
    forager.set_strategies(
            params[Forager::Strategy::s_sigma_A],
            params[Forager::Strategy::s_mean_column_velocity],
            params[Forager::Strategy::s_inspection_time],
            params[Forager::Strategy::s_discrimination_threshold],
            params[Forager::Strategy::s_search_image]
    );
    //ExecutionTimer<std::chrono::milliseconds> nrei_timer("NREI time");
    fitness = forager.NREI();
    //int nrei_duration_ms = nrei_timer.quiet_stop();
    //printf("NREI was %.5f for wolf with params %.5f, %.5f, %.5f, %.5f, %.5f.     TOOK %d ms.\n", double(fitness), params[0], params[1], params[2], params[3], params[4], nrei_duration_ms);
    computing_fitness = false;
}

void Wolf::thread_main() {
    while (thread_running) {
        if (needs_fitness_recalculated) {
            calculate_fitness();
            needs_fitness_recalculated = false;
        }
    }
}

void Wolf::start_thread(){
    thread_running = true;
    the_thread = std::thread(&Wolf::thread_main, this);
}

void Wolf::set_params(Eigen::ArrayXd params) {
    this->params = std::move(params);
    needs_fitness_recalculated = true;
}