//
// Created by Jason Neuswanger on 5/14/18.
//

#ifndef VALIDFISH_WOLF_H
#define VALIDFISH_WOLF_H

#include "Optimizer.h"
#include "Forager.h"

class Optimizer;
class Forager;

class Wolf {

private:

    std::thread the_thread;
    std::atomic<bool> thread_running = false;
    std::atomic<bool> needs_fitness_recalculated = true;  // Initialized to true for initial fitness calculation
    Forager forager;
    Optimizer *optimizer;
    void enforce_bounds_and_context();
    void calculate_fitness();
    void thread_main();

public:

    Eigen::ArrayXd params;
    std::atomic<double> fitness;
    std::atomic<bool> computing_fitness = false;

    bool operator > (const Wolf &other_wolf) const { return (fitness > other_wolf.fitness); }

    Wolf(Optimizer *optimizer, Forager *forager, size_t n_vars);

    ~Wolf();
    void start_thread();
    void set_params(Eigen::ArrayXd params);

};


#endif //VALIDFISH_WOLF_H
