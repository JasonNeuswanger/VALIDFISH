//
// Created by Jason Neuswanger on 2/16/18.
//

#ifndef DRIFTMODELC_SWIMMER_H
#define DRIFTMODELC_SWIMMER_H

#include <algorithm>
#include <cassert>
#include <experimental/filesystem>
#include <map>
#include <memory>

#include "ManeuverInterpolation.h"

namespace fs = std::experimental::filesystem;

class Swimmer {

private:

    double brett_glass_regression_value(unsigned which_params);

    double SMR;                 // Standard metabolic rate (mg O2 * kg/h)
    double AMR;                 // Active metabolic rate (mg O2 * kg/h)
    double u_ms;                // Maximum sustainable swimming speed in m/s

    std::map<double, std::shared_ptr<ManeuverInterpolation>> interps;

protected:  // Functions and attributes available to subclasses of Swimmer, but not outside classes

    void process_parameter_updates(bool rebuild_interpolations);

    fs::path maneuver_interpolation_csv_base_path;
    double mass_g;              // Mass in grams
    double fork_length_cm;      // Fork length in cm
    unsigned temperature_C;     // Temperature in degrees C, rounded to the nearest integer

public:

    Swimmer(double fork_length_cm,
            double mass_g,
            unsigned temperature_C,
            std::string *maneuver_interpolation_csv_base_path);

    Swimmer (Swimmer *otherSwimmer);

    /* After creating a swimmer object, you can't directly modify its mass, temperature, etc., because modifying
     * any of those requires recalculating a bunch of other things. Instead, there's a general setter function
     * set_parameters() to modify any/all parameters at once and then do all the required updates. */

    void modify_parameters(double fork_length_cm,
                           double mass_g,
                           unsigned temperature_C);

    /* steady_swimming_cost() returns the energy cost in J/s of steady swimming at swimming_speed m/s
     * It is calculated from the same Brett & Glass (1973) sockeye salmon model we use for maneuver costs. */

    double steady_swimming_cost(double swimming_speed);

    /* maneuver_cost() returns the interpolated cost (energy in J/s if is_energy_cost==True, or pursuit duration in s
     * if is_energy_cost==False) of maneuvering to capture a prey at (x,y,z) relative to the fish's position (0,0,0).
     * Units for the input coordinates are meters (later converted to cm internally), m/s for velocity (later converted
     * to cm/s internally), and the input coordinates are:
     * x: positive to the right, negative to the left of an upstream-facing fish
     * y: positive upstream, negative downstream
     * z: vertical, positive up, negative down */

    double maneuver_cost(double x, double y, double z, double v, bool is_energy_cost);

};


#endif //DRIFTMODELC_SWIMMER_H
