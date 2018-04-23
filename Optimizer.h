//
// Created by Jason Neuswanger on 2/12/18.
//

#ifndef DRIFTMODELC_OPTIMIZER_H
#define DRIFTMODELC_OPTIMIZER_H

#include "Forager.h"
#include "utility.h"
#include <gsl/gsl_randist.h>

using namespace Eigen;

struct wolf_type {
    ArrayXd params;
    double fitness;
    bool operator > (const wolf_type& other_wolf) const { return (fitness > other_wolf.fitness); }
};

class Optimizer {

private:

    Forager *initial_forager;

    bool verbose;

    ArrayXd r1, r2, A1, A2, A3, C1, C2, C3, D_alpha, D_beta, D_delta, X_1, X_2, X_3;

    std::vector<wolf_type> wolves;
    wolf_type *alpha, *beta, *delta;

    size_t max_iterations, pack_size;

    size_t n_vars;

    std::map<Forager::Strategy, double> context;    // holds any values we want fixed

    ArrayXd random_weights();
    ArrayXd new_weights(ArrayXd previous_weights);
    double validated_random_parameter_value(Forager::Strategy s);
    wolf_type random_wolf();
    double wolf_fitness(wolf_type *wolf);
    void calculate_wolf_fitnesses();
    void update_forager_from_wolf(Forager *forager, wolf_type *wolf);
    void enforce_bounds_and_constraints();


public:

    bool algorithm_use_chaos = false;               // CGWO from Kohli & Arora 2017 (Chaotic GWO)
    bool algorithm_use_dynamic_C = false;           // MGWO from Kumar et al 2016
    bool algorithm_use_exponential_decay = false;   // mGWO from Mittal et al 2016
    bool algorithm_use_levy = false;                // LGWO from Luo et al 2017 (Levy flight GWO)
    bool algorithm_use_only_alpha = false;          // EGWO from Joshi & Arora 2017
    bool algorithm_use_weighted_alpha = true;       // Custom modification of EGWO

    void set_algorithm_options(bool use_chaos, bool use_dynamic_C, bool use_exponential_decay, bool use_levy,
                                   bool use_only_alpha, bool use_weighted_alpha);

    void add_context(Forager::Strategy s, double value);
    void clear_context();

    // one more method I haven't implemented yet: http://journals.sagepub.com/doi/pdf/10.1177/1176934317729413

    Optimizer(Forager *initial_forager, size_t _max_evaluations, size_t pack_size, bool verbose);

    std::vector<double> optimize_forager();

};


#endif //DRIFTMODELC_OPTIMIZER_H
