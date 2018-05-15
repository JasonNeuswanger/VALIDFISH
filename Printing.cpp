//
// Created by Jason Neuswanger on 2/28/18.
//

#include "Forager.h"


void Forager::print_cache_sizes() {
    printf("Memoization cache activity:\n");
    printf("                  Expected maneuver cache hits %10ld      misses %10ld.\n", expected_maneuver_cost_cache_hits, expected_maneuver_cost_cache_misses);
    printf("              Detection probability cache hits %10ld      misses %10ld.\n", detection_probability_cache_hits, detection_probability_cache_misses);
    printf("                Mean value function cache hits %10ld      misses %10ld.\n", mean_value_function_cache_hits, mean_value_function_cache_misses);
    printf("                                Tau cache hits %10ld      misses %10ld.\n", tau_cache_hits, tau_cache_misses);
    printf("                      Detection PDF cache hits %10ld      misses %10ld.\n", detection_pdf_cache_hits, detection_pdf_cache_misses);
    printf("Expected discrimination probability cache hits %10ld      misses %10ld.\n", expected_discrimination_probability_cache_hits, expected_discrimination_probability_cache_misses);
    printf("         Discrimination probability cache hits %10ld      misses %10ld.\n", discrimination_probability_cache_hits, discrimination_probability_cache_misses);
    printf("            Bounds of profitability cache hits %10ld      misses %10ld.\n", bounds_of_profitability_cache_hits, bounds_of_profitability_cache_misses);
    printf("There were %ld numerator evaluations and %ld denominator evaluations.\n",  numerator_integrand_evaluations, denominator_integrand_evaluations);
}

void Forager::print_strategy() {
    printf("Foraging strategy:\n");
    printf("                 sigma_A : %.5f radians\n", sigma_A);
    printf("    mean column velocity : %.5f m/s        (focal velocity %.5f m/s)\n", mean_column_velocity, focal_velocity);
    printf("         inspection time : %.8f s\n", inspection_time);
    printf("discrimination threshold : %.8f\n", discrimination_threshold);
    std::shared_ptr<PreyType> search_image_type = nullptr;
    for (auto & pt : prey_types) {
        if (pt->search_image_status == PreyType::SearchImageStatus::search_image_target) {
            search_image_type = pt;
        }
    }
    if (search_image_type == nullptr) {
        printf("        search image for : None\n");
    } else {
        printf("        search image for : %s\n", search_image_type->name.c_str());
    }
}

void Forager::print_parameters() {
    printf("Parameters:\n");
    for (int pInt = first_parameter; pInt <= last_parameter; pInt++) {
        auto p = static_cast<Parameter>(pInt);
        double value = get_parameter(p);
        if (value > 0.0001 && value < 10000) {
            printf("%20s : %.5f\n", parameter_names[p].c_str(), value);
        } else {
            printf("%20s : %.3e\n", parameter_names[p].c_str(), value);
        }
    }
}

//void Forager::print_discrimination_probabilities() {
//    for (auto & pt : prey_types) {
//        printf("For category %20.20s, p(false_positive)=blank and p(true_hit)=blank. Perceptual sigma=blank (NOT IMPLEMENTED.)\n", pt->name.c_str()); // todo implement discrimination prob printing
//    }
//}

void Forager::print_analytics(){
    analyze_results();
    printf("Overall, pursued %.5f items/s (%.5f prey, %.5f debris). Ingested %.5f of items pursued.\n", foraging_attempt_rate, prey_pursuit_rate, debris_pursuit_rate, proportion_of_attempts_ingested);
    for (auto & pt : prey_types) {
        printf("For %20.20s, pursued %.5f items/s (%.5f prey, %.5f debris). Ingested %.5f of items pursued.\n", pt->name.c_str(), pt->foraging_attempt_rate, pt->prey_pursuit_rate, pt->debris_pursuit_rate, pt->proportion_of_attempts_ingested);
    }
}

void Forager::time_NREIs(size_t iters, size_t nreis_per_iter) {
    double result;
    for (int i=0; i < iters; i++) {
        ExecutionTimer<std::chrono::milliseconds> nrei_timer("ms");
        for (int j=0; j < nreis_per_iter; j++) {
            process_parameter_updates();
            result = NREI();
        }
        printf("Computed %lu NREIs. Latest value %.5f. ", nreis_per_iter, result);
        nrei_timer.stop();
    }
}

void Forager::print_status() {
    double cs_area = cross_sectional_area();
    double max_volume = volume_within_radius(max_radius);
    printf("Search volume for largest prey type is % .5f\n", max_volume);
    printf("Focal velocity is % .5f\n", focal_velocity);
    printf("Cross-sectional area is % .5f\n", cs_area);
    compute_set_size(true); // prints set size substats
    printf("Set size is % .5f\n", set_size);
    printf("Focal swimming cost is % .5f\n", focal_swimming_cost);
    // print_discrimination_probabilities();
    print_analytics();

    // The ones below work okay for the example but get into trouble on real fish data because they go out of bounds.
    // PreyType pc1 = get_prey_type(&pc1_name);
    //    double test_x = 0.05;
//    double test_y = 0.1;
//    double test_z = 0.08;
//    double test_v = 0.20;
//    double detprob = detection_probability(test_x, test_z, &pc1);
//    printf("Detection probability for %s at (% .2f, % .2f) is % .18f\n", pc1.name.c_str(), test_x, test_x, detprob);
//    double energy_cost = maneuver_cost(test_x, test_y, test_z, test_v, true);
//    double pursuit_duration = maneuver_cost(test_x, test_y, test_z, test_v, false);
//    printf("Energy cost at (x=%.2f, y=%.2f, z=%.2f) interpolated to be %.8f J.\n", test_x, test_y, test_z, energy_cost);
//    printf("Pursuit duration at (x=%.2f, y=%.2f, z=%.2f) interpolated to be %.8f s.\n", test_x, test_y, test_z, pursuit_duration);
//    double mmec = expected_maneuver_cost(test_x, test_z, &pc1, true, detprob);
//    double mmpd = expected_maneuver_cost(test_x, test_z, &pc1, false, detprob);
//    printf("Mean maneuver energy cost for type %s at (x=%.2f, z=%.2f) is %.8f J.\n", pc1.name.c_str(), test_x, test_z, mmec);
//    printf("Mean maneuver pursuit duration for type %s at (x=%.2f, z=%.2f) is %.8f J.\n", pc1.name.c_str(), test_x, test_z, mmpd);
    double nrei = NREI();
    printf("NREI is %.8f.\n", nrei);

    /*
    double nrei;
    ExecutionTimer<std::chrono::milliseconds> timer1("10 NREI calculations");
    for (int i = 0; i < 10; ++i) {
        //printf("i=%d\n", i);
        expected_maneuver_cost_cache.clear();   // the clears don't take up much time at all within the loop
        detection_probability_cache.clear();
        nrei = NREI();
    }
    timer1.stop();

    /*

    // SANDBOX FOR SMALL FUNCTION TESTS
    DIAG_NOCACHE = true;
    double mmc, pd, x, z, v, *xz;
    ExecutionTimer<std::chrono::milliseconds> timer1("cost without cache");
    for (int i = 0; i < 1000000; i++) {
        xz = random_xz(); x = xz[0]; z = xz[1];
        //pd = detection_probability(x, x, &pc1);
        //mmc = mean_maneuver_cost(xz[0], xz[1], &pc1, true, pd);
        v = water_velocity(z);
    }
    timer1.stop();
    DIAG_NOCACHE = false;
    ExecutionTimer<std::chrono::milliseconds> timer2("cached cost");
    for (int i = 0; i < 1000000; i++) {
        xz = random_xz(); x = xz[0]; z = xz[1];
        //pd = detection_probability(xz[0], xz[1], &pc1);
        //mmc = expected_maneuver_cost(xz[0], xz[1], &pc1, true, pd);
        v = water_velocity(z);
    }
    timer2.stop();
     */




    //print_cache_sizes();

}
