//
// Created by Jason Neuswanger on 5/13/18.
//

#ifndef VALIDFISH_OPTIMIZER_H
#define VALIDFISH_OPTIMIZER_H

#include "Forager.h"
#include "Wolf.h"
#include "utility.h"
#include <gsl/gsl_randist.h>

using namespace Eigen;
class Wolf;

class Optimizer {

private:

    std::default_random_engine generator;
    std::unordered_map<Forager::Strategy, std::uniform_real_distribution<double>> distributions;

    bool verbose;

    ArrayXd r1, r2, A1, A2, A3, C1, C2, C3, D_alpha, D_beta, D_delta, X_1, X_2, X_3;

    std::vector<std::shared_ptr<Wolf>> wolves;
    Wolf *alpha, *beta, *delta;

    size_t max_iterations;
    size_t n_vars;

    std::unordered_map<Forager::Strategy, double> context;  // Holds values we want to fix instead of optimize

    ArrayXd random_weights();

    void wait_for_wolves_to_calculate_fitnesses();


public:

    std::shared_ptr<Forager> initial_forager;
    double validated_random_strategy_value(Forager::Strategy s);

    void add_context(Forager::Strategy s, double value);
    double context_value(Forager::Strategy s);
    void clear_context();

    Optimizer(std::shared_ptr<Forager> initial_forager, size_t max_iterations, size_t pack_size, bool verbose);

    std::vector<double> optimize_forager();


};


#endif //VALIDFISH_OPTIMIZER_H
