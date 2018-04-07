//
// Created by Jason Neuswanger on 1/31/18.
//

#ifndef DRIFTMODELC_FORAGER_H
#define DRIFTMODELC_FORAGER_H

#include <functional>
#include <vector>
#include <cassert>
#include <iostream>
#include <experimental/filesystem>
#include <thread>
#include <future>

#include <Eigen/Dense>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_cdf.h>

#include "third-party/gsl_function_pp.h"
#include "third-party/ExecutionTimer.h"
#include "third-party/flat_hash_map.hpp"

#include "Swimmer.h"
#include "PreyCategory.h"
#include "utility.h"

#define USE_ADAPTIVE_INTEGRATION false  // Set to false for the fastest algorithm and general use
#define QUAD_SUBINT_LIM 100             // only relevant when using adaptive integration
#define QUAD_EPSABS 0                   // always set to 0 to use relative instead
#define QUAD_EPSREL 1e-0                // set precision, from 1e0 to 1e-3 or so (fast to slow)
#define MEMOIZATION_PRECISION 1e-2      // spatial precision of the memoized caches, should be 0.01*QUAD_EPSREL or better
#define GSL_ERROR_POLICY 0              // 0 to ignore GSL errors. 1 to print them but continue. 2 to abort (for debugging mode).

class Forager : public Swimmer {

private:

    // Diagnostics / performance tracking

    bool DIAG_NOCACHE = false;
    bool DIAG_NANCHECKS = true;
    size_t denominator_integrand_evaluations = 0;
    size_t numerator_integrand_evaluations = 0;

    // Fish decision variables specified as the main part of its strategy (other than prey attention allocations)
    //double radius, theta, mean_column_velocity, saccade_time, discrimination_threshold;

    // Quantities derived from the above but saved as instance variables rather than recalculated whenever needed
    double search_volume, set_size, search_rate, focal_velocity, focal_swimming_cost, angular_resolution;

    // Variables pertaining to characteristics of the environment
    double surface_z, bottom_z, depth, bed_roughness, lambda_c, sigma_t, base_crypticity;

    // Model parameters
    double delta_0; // Scales the effect of angular size on tau
    double alpha_0; // Scales the effect of feature-based attention on tau
    double beta;    // Scales the effect of saccade time * set size on tau
    double Z_0;     // Scales the effect of spatial attention on tau; viewable as attentional capacity
    double c_1;     // Scales effect of feature-based attention and saccade time on discrimination

    // Model design choices
    double exclude_unprofitable_maneuvers = true;

    // Functions in Forager.cpp

    void process_parameter_updates();       // Handles any change to parameters or strategy variables -- MAKE PRIVATE EVENTUALLY
    double RateOfEnergyIntake(bool is_net); // Pass is_net = true for NREI, false for GREI
    void set_bounds();                      // Initialization step for optimization bounds
    void print_cache_sizes();               // Only used for diagnostics

    // Geometry.cpp

    double integrate_over_xz_plane(gsl_function *, bool integrand_is_1d);
    double integrate_over_volume(gsl_function *func, double min_rho, double max_rho, double min_theta, double max_theta);
    double* random_xz();

    // EnvironmentInhabitant.cpp

    std::vector<PreyCategory> prey_categories;

    void normalize_feature_sizes();

    // EnergyCosts.cpp

    double integrate_detection_pdf(double x, double z, PreyCategory *pc);
    double integrate_energy_cost_over_prey_path(double x, double z, PreyCategory *pc, bool is_energy_cost);

    // General recalculating

    void compute_angular_resolution();      // Forager.cpp
    void compute_set_size(bool verbose);    // Forager.cpp
    void compute_search_rate();             // Forager.cpp
    void compute_focal_velocity();          // EnvironmentInhabitant.cpp
    void compute_discrimination_probabilities();    // Forager.cpp

    // Caches to save the results of expensive calculations that are repeated exactly (i.e., memoization)

    template< typename hash_key_type >
    struct hash_function {
        hash_key_type operator()(const hash_key_type &key) {
            return key;
        }
    };
    ska::flat_hash_map<long long, double, hash_function<long long>> mean_maneuver_cost_cache;
    ska::flat_hash_map<long long, double, hash_function<long long>> detection_probability_cache;
    size_t mean_maneuver_cost_cache_hits = 0;
    size_t mean_maneuver_cost_cache_misses = 0;
    size_t detection_probability_cache_hits = 0;
    size_t detection_probability_cache_misses = 0;

    // Analysis/saving of results

    double prey_pursuit_rate, debris_pursuit_rate, foraging_attempt_rate, proportion_of_attempts_ingested;

    double pursuit_rate(std::string which_rate, PreyCategory *pc);

    // Optimization

    typedef std::unordered_map<std::string, double[2]> BoundsMap;

public:

    BoundsMap bounds;

    double radius, theta, mean_column_velocity, saccade_time, discrimination_threshold; // TEMPORARILY PUBLIC

    // Forager.cpp

    Forager(double fork_length_cm,
            double mass_g,
            double radius,
            double theta,
            double mean_column_velocity,
            double saccade_time,
            double discrimination_threshold,
            double delta_0,
            double alpha_0,
            double Z_0,
            double c_1,
            double beta,
            double bottom_z,
            double surface_z,
            unsigned temperature,
            double bed_roughness,
            double lambda_c,
            double sigma_t,
            double base_crypticity,
            std::string *maneuver_interpolation_csv_base_path);
    Forager(Forager *otherForager);     // Copy constructor
    virtual ~Forager();                 // Destructor

    void modify_strategies(double radius, double theta, double mean_column_velocity, double saccade_time, double discrimination_threshold, std::vector<double> attention);
    void modify_parameters(double delta_0, double alpha_0, double beta, double Z_0, double c_1);
    enum Strategy { s_radius, s_theta, s_mean_column_velocity, s_saccade_time, s_discrimination_threshold };
    enum Parameter { p_delta_0, p_alpha_0, p_beta, p_Z_0, p_c_1 };
    void modify_strategy(Strategy strategy, double value);
    void modify_parameter(Parameter parameter, double value);
    void modify_bound(std::string field, double lower_bound, double upper_bound);
    void fix_bound(std::string field, double fixed_value);

    double tau(double t, double x, double z, PreyCategory *pc);
    double detection_probability(double x, double z, PreyCategory *pc);
    double mean_maneuver_cost(double x, double z, PreyCategory *pc, bool is_energy_cost, double det_prob);

    double NREI();    // Net rate of energy intake
    double GREI();    // Gross rate of energy intake

    // Printing.cpp

    void print_status();
    void print_strategy();
    void print_discrimination_probabilities();
    void print_analytics();

    // ForagerGetters.cpp

    PreyCategory *get_prey_category(std::string *name);
    double get_fork_length_cm();
    double get_radius();
    double get_theta();
    double get_focal_velocity();
    double get_foraging_attempt_rate();
    double get_proportion_of_attempts_ingested();
    double get_diet_proportion_for_prey_category(std::string *category_name);

    // Geometry.cpp

    double cross_sectional_area();
    double volume_within_radius(double r);

    // EnvironmentInhabitant.cpp

    double water_velocity(double z);
    void evenly_distribute_attention();
    void build_sample_prey_categories();
    size_t num_prey_categories();
    void add_prey_category(int number, std::string name, double mean_prey_length, double mean_prey_energy,
                                    double crypticity_multiplier, double prey_drift_density, double debris_drift_density, double feature_size);
    void process_prey_category_changes();

    // Analysis.cpp -- functions for describing the forager for model testing.

    // These generally take prey categories by name instead of the category object, which makes it easier to call them
    // for a given prey category in Python without having to instantiate the prey category objects in Python or worry
    // about whether they're pointing to the same object or a copy.

    void analyze_results();
    double depletable_detection_probability(double x, double y, double z, std::string *pc_name);
    double relative_pursuits_by_position_single_prey_type(double x, double y, double z, std::string *pc_name);
    double relative_pursuits_by_position(double x, double y, double z); // sums the above over all prey categories
    double proportion_of_detections_within(double min_distance, double max_distance, double min_angle, double max_angle,
                                               std::string *pc_name, std::string *which_items);

};

#endif //DRIFTMODELC_FORAGER_H
