//
// Created by Jason Neuswanger on 2/28/18.
//

#include "Forager.h"


void Forager::print_cache_sizes() {
    printf("Mean maneuver cache hits %ld misses %ld. Detection probability hits %ld misses %ld.\n",
           mean_maneuver_cost_cache_hits, mean_maneuver_cost_cache_misses,
           detection_probability_cache_hits, detection_probability_cache_misses);
    printf("There were %ld numerator evaluations and %ld denominator evaluations.\n", numerator_integrand_evaluations, denominator_integrand_evaluations);
}

void Forager::print_strategy() {
    printf("Foraging strategy:\n");
    printf("radius                   : %.5f m\n", radius);
    printf("theta                    : %.5f radians\n", theta);
    printf("mean column velocity     : %.5f m/s        (focal velocity %.5f m/s)\n", mean_column_velocity, focal_velocity);
    printf("saccade time             : %.8f s\n", saccade_time);
    printf("discrimination threshold : %.8f\n", discrimination_threshold);
    for (auto &pc : prey_categories) {
        printf("attention to category %s : %.5f\n", pc.name.c_str(), pc.get_attention_allocated());
    }
}

void Forager::print_discrimination_probabilities() {
    for (auto &pc : prey_categories) {
        printf("For category %20.20s, p(false_positive)=%.8f and p(true_hit)=%.8f.\n", pc.name.c_str(), pc.false_positive_probability, pc.true_hit_probability);
    }
}

void Forager::print_analytics(){
    analyze_results();
    printf("Overall, pursued %.5f items/s (%.5f prey, %.5f debris). Ingested %.5f of items pursued.\n", foraging_attempt_rate, prey_pursuit_rate, debris_pursuit_rate, proportion_of_attempts_ingested);
    for (auto &pc : prey_categories) {
        printf("For %20.20s, pursued %.5f items/s (%.5f prey, %.5f debris). Ingested %.5f of items pursued.\n", pc.name.c_str(), pc.foraging_attempt_rate, pc.prey_pursuit_rate, pc.debris_pursuit_rate, pc.proportion_of_attempts_ingested);
    }
}

void Forager::print_status() {
    std::string pc1_name = "4-5 mm size class";
    double cs = cross_sectional_area();
    printf("Search volume is % .18f\n", search_volume);
    printf("Focal velocity is % .18f\n", focal_velocity);
    printf("Cross-sectional area is % .18f\n", cs);
    printf("Search rate is % .18f\n", search_rate);
    compute_set_size(true); // prints set size substats
    printf("Set size is % .18f\n", set_size);
    printf("Focal swimming cost is % .18f\n", focal_swimming_cost);
    print_discrimination_probabilities();
    print_analytics();

    // The ones below work okay for the example but get into trouble on real fish data because they go out of bounds.
    // PreyCategory pc1 = get_prey_category(&pc1_name);
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
//    double mmec = mean_maneuver_cost(test_x, test_z, &pc1, true, detprob);
//    double mmpd = mean_maneuver_cost(test_x, test_z, &pc1, false, detprob);
//    printf("Mean maneuver energy cost for type %s at (x=%.2f, z=%.2f) is %.8f J.\n", pc1.name.c_str(), test_x, test_z, mmec);
//    printf("Mean maneuver pursuit duration for type %s at (x=%.2f, z=%.2f) is %.8f J.\n", pc1.name.c_str(), test_x, test_z, mmpd);
    double nrei = NREI();
    printf("NREI is %.8f.\n", nrei);

    /*
    double nrei;
    ExecutionTimer<std::chrono::milliseconds> timer1("10 NREI calculations");
    for (int i = 0; i < 10; ++i) {
        //printf("i=%d\n", i);
        mean_maneuver_cost_cache.clear();   // the clears don't take up much time at all within the loop
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
        //mmc = mean_maneuver_cost(xz[0], xz[1], &pc1, true, pd);
        v = water_velocity(z);
    }
    timer2.stop();
     */




    //print_cache_sizes();

}
