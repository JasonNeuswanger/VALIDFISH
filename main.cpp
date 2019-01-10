#include <iostream>

#include "Forager.h"
#include "Optimizer.h"
#include "Wolf.h"

int main() {

    std::string maneuver_interpolation_path = "/Users/Jason/Dropbox/Drift Model Project/Calculations/driftmodeldev/maneuver-model-tables/";

    // Sample fish based on 2015-07-17-3 Panguingue - Dolly Varden (id #4)
    auto forager = std::make_shared<Forager>(  14.34,     // double fork_length_cm,
                            23.28,  // double mass_g,
                            1.0,    // double sigma_A
                            0.59,   // double mean_column_velocity,
                            0.5,    // double inspection_time,
                            2.0,    // double discrimination_threshold,
                            -1,     // double search_image,
                            1e-4,   // double delta_0,
                            10,     // double alpha_tau,
                            10,     // double alpha_d,
                            0.5,    // double A_0,
                            0.2,    // double beta,
                            -0.04,  // double bottom_z,
                            0.47,   // double surface_z,
                            10,     // unsigned temperature,
                            0.5,    // double tau_0,
                            50,     // double flicker_frequency,
                            1e-3,   // double nu_0,
                            2.0,    // double discriminability,
                            1e-2,   // double delta_p,
                            1.0,    // double omega_p,
                            0.5,    // double ti_p,
                            1.0,    // double sigma_p_0
                            &maneuver_interpolation_path);
    forager->build_sample_prey_types();

//    ExecutionTimer<std::chrono::milliseconds> nrei_timer("Single NREI time");
//    printf("Forager NREI is %.6f.\n", forager->NREI());
//    nrei_timer.stop();
//    forager->time_NREIs(5, 3);
//    forager->print_cache_sizes();

    // auto test_pt = forager->get_prey_type("2 mm class");

    //forager->spatial_detection_proportions(test_pt, "All", true);

//    ExecutionTimer<std::chrono::milliseconds> opt_timer("Optimization time: ");
//    Optimizer *optimizer;
//    bool verbose = true;
//    size_t n_iterations = 20;
//    size_t pack_size = 6;
//    optimizer = new Optimizer(forager, n_iterations, pack_size, verbose);
//    optimizer->optimize_forager();
//    opt_timer.stop();

    // Optimization times vs pack sizes (20 iterations):
    // Pack size 6: 22313, 24405, 24923
    // Pack size 7: 38070, 26931, 27718
    // Pack size 8: 35740, 26585, 40418
    // Pack size 9: 31502, 34441, 33953
    // It looks like execution times vary quite a bit depending on what other threads on the system are doing,
    // but the most reliable performance comes from using 2 less than the number of processor cores.

    forager->print_prey();
    forager->print_strategy();
    forager->print_parameters();
    forager->print_analytics();
    printf("\nForager NREI is %.6f.\n", forager->NREI());

}
