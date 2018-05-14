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
#include <gsl/gsl_randist.h>

#include "lib/gsl_function_pp.h"
#include "lib/ExecutionTimer.h"
#include "lib/flat_hash_map.hpp"

#include "Swimmer.h"
#include "PreyType.h"
#include "utility.h"

#define USE_ADAPTIVE_INTEGRATION false  // Set to false for the fastest algorithm and general use
#define QUAD_SUBINT_LIM 100             // only relevant when using adaptive integration
#define QUAD_EPSABS 0                   // always set to 0 to use relative instead
#define QUAD_EPSREL 1e0                 // set precision, from 1e0 to 1e-3 or so (fast to slow) -- 1e0 works fine
#define MEMOIZATION_PRECISION 0.015     // spatial precision of the memoized caches in m, currently 1.5 cm
#define MEMOIZATION_PRECISION_TIME 0.1  // temporal precision of memoized caches, in seconds
#define MEMOIZATION_PRECISION_PROBABILITY 0.01  // temporal precision of detection probability based cache
#define GSL_ERROR_POLICY 0              // 0 to ignore GSL errors. 1 to print them but continue. 2 to abort (for debugging mode).

#define DIAG_NOCACHE false              // Disable various internal caches (works MUCH more slowly if set to true)
#define DIAG_NOCACHE_TAU false
#define DIAG_NOCACHE_DETECTION_PROBABILITY false
#define DIAG_NOCACHE_DETECTION_PDF false
#define DIAG_NOCACHE_MEAN_VALUE_FUNCTION false   // Given the other caches, this one was just a pointless lookup that almost never hit

#define DIAG_NO_DISCRIMINATION_MODEL false  // Discrimination model takes up about half the time; that's where the savings are

class Forager : public Swimmer {

private:

    // Performance tracking
    size_t denominator_integrand_evaluations = 0;
    size_t numerator_integrand_evaluations = 0;

    // Quantities derived from the above but saved as instance variables rather than recalculated whenever needed
    double set_size, focal_velocity, focal_swimming_cost, angular_resolution;

    // Variables pertaining to characteristics of the environment
    double surface_z, bottom_z, depth, discriminability, flicker_frequency;

    // Model parameters for detection
    double tau_0;                   // Bare minimum value all the other tau factors multiply
    double delta_0;                 // Scales the effect of angular size on tau
    double alpha_tau;               // Scales the effect of search image, if any, on tau
    double beta;                    // Scales the effect of inspection time * set size on tau
    double A_0;                     // Scales the effect of spatial attention on tau; viewable as attentional capacity
    double nu_0;                    // Scales effect of loom on tau
    // Model parameters for discrimination
    double delta_p;                 // Scales effect of angular size on perceptual variance
    double omega_p;                 // Scales effect of angular velocity on perceptual variance
    double ti_p;                    // Scales effect of inspection time on perceptual variance
    double sigma_p_0;               // Base level of perceptual variance, whihc all other factors multiply.
    double alpha_d;                 // Scales the effect of search image, if any, on discrimination

    // Scales relative effect of existing search volume vs incoming water volume on search rate

    // Physiological constraints borrowed from the literature
    double min_prey_length_from_gill_rakers, max_prey_length_from_mouth_gape;
    double max_radius = 0;  // shorthand for the visible radius of the largest prey type

    // Functions in Forager.cpp

    void process_parameter_updates();       // Handles any change to parameters or strategy variables -- MAKE PRIVATE EVENTUALLY
    double RateOfEnergyIntake(bool is_net, bool is_cost); // Pass is_net = true for NREI, false for GREI
    void set_parameter_bounds();            // Ditto for parameter optimization bounds
    void set_strategy_bounds();

    // DetectionSubmodel.cpp

    void compute_set_size(bool verbose);    // Forager.cpp

    // Private internal methods here calculate the quantities below, whereas the public methods check the cache and call these if needed.
    double calculate_mean_value_function(double T, double x, double z, const PreyType &pt);
    double calculate_detection_probability(double x, double z, const PreyType &pt);


    // DiscriminationSubmodel.cpp

    double average_discrimination_probability_over_prey_path(double x, double z, const PreyType &pt,
                                                             bool is_false_positive, double det_prob);
    std::pair<double, double> discrimination_probabilities(double t, double x, double z, const PreyType &pt);

    std::pair<double, double> expected_discrimination_probabilities(double x, double z, const PreyType &pt,
                                                                    double det_prob);
    double perceptual_variance(double t, double x, double z, const PreyType &pt);
    inline double perception_effect_of_angular_area(double distance, const PreyType &pt);
    inline double perception_effect_of_angular_velocity(double v, double t, double xsq, double zsq, double rsq);
    inline double perception_effect_of_search_image(const PreyType &pt);
    inline double perception_effect_of_inspection_time();

    // Geometry.cpp

    double integrate_over_xz_plane(gsl_function *func, bool integrand_is_1d);
    double integrate_over_volume(gsl_function *func, double min_rho, double max_rho, double min_theta, double max_theta);
    double time_at_y(double y, double x, double z, const PreyType &pt);
    double y_at_time(double t, double x, double z, const PreyType &pt);
    double* random_xz();

    // EnvironmentInhabitant.cpp

    std::vector<std::shared_ptr<PreyType>> prey_types;

    // Energetics.cpp

    double integrate_energy_cost_over_prey_path(double x, double z, const PreyType &pt, bool is_energy_cost);

    // General recalculating

    void compute_angular_resolution();      // Forager.cpp
    void compute_focal_velocity();          // EnvironmentInhabitant.cpp

    // Caches to save the results of expensive calculations that are repeated exactly (i.e., memoization)

    template< typename hash_key_type >
    struct hash_function {
        hash_key_type operator()(const hash_key_type &key) {
            return key;
        }
    };
    ska::flat_hash_map<long long, double, hash_function<long long>> expected_maneuver_cost_cache;
    ska::flat_hash_map<long long, double, hash_function<long long>> detection_probability_cache;
    ska::flat_hash_map<long long, double, hash_function<long long>> mean_value_function_cache;
    ska::flat_hash_map<long long, double, hash_function<long long>> tau_cache;
    ska::flat_hash_map<long long, double, hash_function<long long>> detection_pdf_cache;
    ska::flat_hash_map<long long, std::pair<double, double>, hash_function<long long>> expected_discrimination_probability_cache;
    ska::flat_hash_map<long long, std::pair<double, double>, hash_function<long long>> discrimination_probability_cache;
    ska::flat_hash_map<long long, std::pair<double, double>, hash_function<long long>> bounds_of_profitability_cache;
    size_t expected_maneuver_cost_cache_hits = 0;
    size_t expected_maneuver_cost_cache_misses = 0;
    size_t detection_probability_cache_hits = 0;
    size_t detection_probability_cache_misses = 0;
    size_t mean_value_function_cache_hits = 0;
    size_t mean_value_function_cache_misses = 0;
    size_t expected_discrimination_probability_cache_hits = 0;
    size_t expected_discrimination_probability_cache_misses = 0;
    size_t discrimination_probability_cache_hits = 0;
    size_t discrimination_probability_cache_misses = 0;
    size_t tau_cache_hits = 0;
    size_t tau_cache_misses = 0;
    size_t detection_pdf_cache_hits = 0;
    size_t detection_pdf_cache_misses = 0;
    size_t bounds_of_profitability_cache_hits = 0;
    size_t bounds_of_profitability_cache_misses = 0;

    // Analysis/saving of results

    double prey_pursuit_rate, debris_pursuit_rate, foraging_attempt_rate, proportion_of_attempts_ingested;

    double pursuit_rate(std::string which_rate, std::shared_ptr<PreyType> pt);

public:

    enum Strategy { s_sigma_A, s_mean_column_velocity, s_inspection_time, s_discrimination_threshold, s_search_image };   // attention included only to specify bounds
    enum Parameter { p_delta_0, p_alpha_tau, p_alpha_d, p_beta, p_A_0, p_discriminability, p_flicker_frequency, p_tau_0, p_nu_0, p_delta_p, p_omega_p, p_ti_p, p_sigma_p_0};

    // Name map so feedback about parameters (validation problems, etc) can print out the actual name and not just a number.
    std::map<Strategy, std::string> strategy_names = {
            {s_sigma_A, "sigma_A"},
            {s_mean_column_velocity, "mean_column_velocity"},
            {s_inspection_time, "inspection_time"},
            {s_discrimination_threshold, "discrimination_threshold"},
            {s_search_image, "search_image"}
    };
    std::map<Parameter, std::string> parameter_names = {
            {p_delta_0, "delta_0"},
            {p_alpha_tau, "alpha_tau"},
            {p_alpha_d, "alpha_d"},
            {p_beta, "beta"},
            {p_A_0, "A_0"},
            {p_discriminability, "discriminability"},
            {p_flicker_frequency, "flicker_frequency"},
            {p_tau_0, "tau_0"},
            {p_nu_0, "nu_0"},
            {p_delta_p, "delta_p"},
            {p_omega_p, "omega_p"},
            {p_ti_p, "ti_p"},
            {p_sigma_p_0, "sigma_p_0"}
    };

    double validate(Strategy s, double v);
    double validate(Parameter p, double v);
    // todo build manipulations
    enum Manipulation { m_prey_multiplier, m_debris_multiplier, m_crypticity_multiplier };
    typedef std::unordered_map<Strategy, std::array<double, 2>> StrategyBoundsMap;
    typedef std::unordered_map<Parameter, std::array<double, 2>> ParameterBoundsMap;
    StrategyBoundsMap strategy_bounds;
    ParameterBoundsMap parameter_bounds;

    double theta = 1.83 * M_PI; // gives a 30 degree blind spot in the rear, Rountree 2009 -- but todo find a better source
    double sigma_A, mean_column_velocity, inspection_time, discrimination_threshold, search_image; // TEMPORARILY PUBLIC, SHOULD MAKE PRIVATE

    // Forager.cpp

    Forager(double fork_length_cm, double mass_g, double sigma_A, double mean_column_velocity,
                double inspection_time, double discrimination_threshold, double search_image, double delta_0, double alpha_tau,
                double alpha_d, double A_0, double beta, double bottom_z, double surface_z,
                unsigned temperature, double tau_0, double flicker_frequency, double nu,
                double discriminability, double delta_p, double omega_p, double ti_p, double sigma_p_0, std::string *maneuver_interpolation_csv_base_path);
    Forager(Forager *otherForager);     // Copy constructor
    virtual ~Forager();                 // Destructor

    void set_strategies(double sigma_A, double mean_column_velocity, double inspection_time,
                        double discrimination_threshold, double search_image);
    void set_parameters(double delta_0, double alpha_tau, double alpha_d, double beta, double A_0,
                        double flicker_frequency, double tau_0, double nu_0,
                        double discriminability, double delta_p, double omega_p, double ti_p, double sigma_p_0);

    void set_strategy(Strategy strategy, double value);
    void set_parameter(Parameter parameter, double value);
    void set_single_strategy_bounds(Strategy strategy, double lower_bound, double upper_bound);
    void fix_single_strategy_bound(Strategy strategy, double fixed_value);

    double NREI();                  // Net rate of energy intake
    double GREI();                  // Gross rate of energy intake
    double maneuver_cost_rate();    // Energy expenditure on maneuvers, per unit time

    // DetectionSubmodel.cpp

    inline double tau_effect_of_spatial_attention(double y, double distance);
    inline double tau_effect_of_set_size();
    inline double tau_effect_of_angular_area(double distance, const PreyType &pt);
    inline double tau_effect_of_loom(double distance, double v, double y, const PreyType &pt);
    inline double tau_effect_of_search_image(const PreyType &pt);

    double tau(double t, double x, double z, const PreyType &pt);
    double calculate_tau(double t, double x, double z, const PreyType &pt);
    double detection_probability(double x, double z, const PreyType &pt);
    double mean_value_function(double T, double x, double z, const PreyType &pt);
    double calculate_detection_pdf_at_t(double t, double x, double z, const PreyType &pt);
    double detection_pdf_at_t(double t, double x, double z, const PreyType &pt);
    double detection_pdf_at_y(double y, double x, double z, const PreyType &pt);

    std::map<std::string, double> tau_components(double t, double x, double z, std::shared_ptr<PreyType> pt); // for diagnostics/plotting

    // Energetics.cpp

    double expected_maneuver_cost(double x, double z, const PreyType &pt, bool is_energy_cost,
                                  double det_prob);
    std::pair<double, double>calculate_bounds_of_profitability(double x, double z, const PreyType &pt);
    inline double item_profitability_at_time(double t, double x, double y, double z, const double maneuver_v, const PreyType &pt);
    std::pair<double, double>bounds_of_profitability(double x, double z, const PreyType &pt);
    bool location_is_profitable(double x, double y, double z, const PreyType &pt);

    // Printing.cpp

    void print_status();
    void print_strategy();
    void print_parameters();
    void print_analytics();
    void print_cache_sizes();
    void time_NREIs(size_t iters, size_t nreis_per_iter);

    // ForagerGetters.cpp

    std::shared_ptr<PreyType> get_prey_type(std::string name);
    std::vector<std::shared_ptr<PreyType>> get_prey_types();
    double get_fork_length_cm();
    double get_max_radius();
    double get_focal_velocity();
    double get_focal_swimming_cost();
    double get_field_of_view();
    double get_foraging_attempt_rate();
    double get_proportion_of_attempts_ingested();
    double get_diet_proportion_for_prey_type(std::shared_ptr<PreyType> pt);
    double get_angular_resolution();
    double get_strategy(Strategy strategy);
    double get_parameter(Parameter parameter);
    std::array<double, 2> get_strategy_bounds(Strategy strategy);
    std::array<double, 2> get_parameter_bounds(Parameter parameter);

    // Geometry.cpp

    double passage_time(double x, double z, const PreyType &pt);
    double cross_sectional_area();
    double volume_within_radius(double r);

    // EnvironmentInhabitant.cpp

    double water_velocity(double z);

    void build_sample_prey_types();
    size_t num_prey_types();
    size_t num_search_image_eligible_prey_types();
    void add_prey_type(int number, std::string name, double mean_prey_length, double mean_prey_energy,
                       double crypticity, double prey_drift_concentration, double debris_drift_concentration,
                       bool search_image_eligible);
    void process_prey_type_changes();

    // Analysis.cpp -- functions for describing the forager for model testing.

    // These generally take prey types by name instead of the type object, which makes it easier to call them
    // for a given prey type in Python without having to instantiate the prey type objects in Python or worry
    // about whether they're pointing to the same object or a copy.

    void analyze_results();
    std::shared_ptr<PreyType> get_favorite_prey_type();
    double relative_pursuits_by_position_single_prey_type(double x, double y, double z, std::shared_ptr<PreyType> pt);
    double relative_pursuits_by_position(double x, double y, double z); // sums the above over all prey categories
    double depleted_prey_concentration_single_prey_type(double x, double y, double z, std::shared_ptr<PreyType> pt);   // items/m3
    double depleted_prey_concentration_total_energy(double x, double y, double z);                                     // J/m3, summing over prey types
    std::map<std::string, std::vector<std::map<std::string, double>>> spatial_detection_proportions(std::shared_ptr<PreyType> pt, std::string which_items, bool verbose);

};

#endif //DRIFTMODELC_FORAGER_H
