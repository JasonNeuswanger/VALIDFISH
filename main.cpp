#include <iostream>

#include "Forager.h"
#include "Optimizer.h"

int main() {

    std::string maneuver_interpolation_path = "/Users/Jason/Dropbox/Drift Model Project/Calculations/driftmodeldev/maneuver-model-tables/";

    Forager *forager;

    forager = new Forager(  -1,     // double fork_length_cm,
                            0.8,    // double mass_g,
                            1.0,    // double sigma_A
                            0.3,   // double mean_column_velocity,
                            0.5,    // double inspection_time,
                            2.0,    // double discrimination_threshold,
                            -1,     // double search_image,
                            1e-4,   // double delta_0,
                            10,     // double alpha_tau,
                            10,     // double alpha_d,
                            0.5,    // double A_0,
                            0.5,    // double beta,
                            -0.1,   // double bottom_z,
                            0.3,    // double surface_z,
                            11,     // unsigned temperature,
                            0.5,    // double tau_0,
                            50,     // double flicker_frequency,
                            1e-3,   // double nu_0,
                            2.0,    // double discriminability,
                            1e-2,   // double delta_p,
                            3.0,    // double omega_p,
                            0.5,    // double ti_p,
                            1.0,    // double sigma_p_0
                            &maneuver_interpolation_path);
    forager->build_sample_prey_types();
    ExecutionTimer<std::chrono::milliseconds> opt_timer("Single NREI time");
    printf("Forager NREI is %.6f.\n", forager->NREI());
    opt_timer.stop();

    // Stop Elapsed: 5336

    forager->print_cache_sizes();

    // NREI fairly recently: 0.026620

    // Result with NO caching whatsoever.
    //    Forager NREI is 0.027032.
    //    Stop Elapsed: 53769 (Single NREI time)


    // auto test_pt = forager->get_prey_type("2 mm class");

    //forager->spatial_detection_proportions(test_pt, "All", true);

//    ExecutionTimer<std::chrono::milliseconds> opt_timer("Optimization time: ");
//    Optimizer *optimizer;
//    bool verbose = true;
//    size_t n_iterations = 50;
//    size_t pack_size = 7;
//    optimizer = new Optimizer(forager, n_iterations, pack_size, verbose);
//    optimizer->set_algorithm_options(false, false, false, false, false, true);
//    optimizer->optimize_forager();
//    opt_timer.stop();

    // todo: wolf fitnesses are all over the map
    // todo: getting nan pursuits

    forager->print_strategy();
    forager->print_status();

}
