#include <iostream>

#include "Forager.h"
#include "Optimizer.h"

int main() {

    std::string maneuver_interpolation_path = "/Users/Jason/Dropbox/Drift Model Project/Calculations/driftmodeldev/maneuver-model-tables/";

    Forager *forager;
    forager = new Forager(-1,     // fork length in cm; set to -1 to use a length-mass regression instead
                    0.8,    // mass
                    0.2,    // radius
                    3.0,    // theta
                    0.2,    // mean_column_velocity
                    0.3,    // saccade_time
                    0.8,    // discrimination_threshold
                    1.0,    // delta_0
                    0.1,    // alpha_0
                    0.003,  // Z_0
                    0.3,    // c_1
                    0.1,    // beta
                   -0.1,    // bottom_z
                    0.3,    // surface_z
                    11,     // temperature (integer)
                    0.05,   // bed_roughness
                    1.0,    // discriminability
                    0.5,    // sigma_t
                    0.3,    // tau_0
                    &maneuver_interpolation_path
    );
    forager->build_sample_prey_categories();
    printf("Forager NREI is %.6f.\n", forager->NREI());

    ExecutionTimer<std::chrono::milliseconds> opt_timer("Optimization time: ");

    Optimizer *optimizer;
    bool verbose = true;
    int n_iterations = 30;
    int pack_size = 15;
    optimizer = new Optimizer(forager, n_iterations, pack_size, verbose);
    optimizer->optimize_forager();
    opt_timer.stop();

    forager->print_strategy();
    forager->print_status();

}
