#include <iostream>

#include "Forager.h"
#include "Optimizer.h"

int main() {

    std::string maneuver_interpolation_path = "/Users/Jason/Dropbox/Drift Model Project/Calculations/driftmodeldev/maneuver-model-tables/";

    Forager *forager;

//    X = delta_0=1.030e-08, beta=1.295e-01, theta_0=4.827, t_s_0=1.061, discriminability=3.044, flicker_frequency=0.100, tau_0=7.787, nu=1.684e-12.


    forager = new Forager(  -1,     // double fork_length_cm,
                            0.8,    // double mass_g,
                            0.02,   // double delta_min,
                            1.0,    // double sigma_A
                            0.20,   // double mean_column_velocity,
                            0.2,    // double saccade_time,
                            2.0,    // double discrimination_threshold,
                            -1,     // double search_image,
                            1e-4,   // double delta_0,
                            10,     // double alpha_tau,
                            10,     // double alpha_d,
                            0.5,    // double A_0,
                            0.1,    // double t_s_0,
                            0.5,    // double beta,
                            -0.1,   // double bottom_z,
                            0.3,    // double surface_z,
                            11,     // unsigned temperature,
                            0.05,   // double bed_roughness,
                            2.0,    // double discriminability,
                            0.1,    // double tau_0,
                            50,     // double flicker_frequency,
                            1e-3,  //double nu,
                            &maneuver_interpolation_path);
    forager->build_sample_prey_types();
    printf("Forager NREI is %.6f.\n", forager->NREI());

    // current test value: Forager NREI is 0.025479.

    auto test_pt = forager->get_prey_type("2 mm class");

    forager->spatial_detection_proportions(test_pt, "All", true);

//    ExecutionTimer<std::chrono::milliseconds> opt_timer("Optimization time: ");
//
//    Optimizer *optimizer;
//    bool verbose = true;
//    size_t n_iterations = 50;
//    size_t pack_size = 7;
//    optimizer = new Optimizer(forager, n_iterations, pack_size, verbose);
//    optimizer->set_algorithm_options(false, false, false, false, false, true);
//    optimizer->optimize_forager();
//    opt_timer.stop();
    forager->print_strategy();
    forager->print_status();

}
