//
// Created by Jason Neuswanger on 2/19/18.
//

#ifndef DRIFTMODELC_MANEUVERINTERPOLATION_H
#define DRIFTMODELC_MANEUVERINTERPOLATION_H

#include <algorithm>
#include <experimental/filesystem>

#include <gsl/gsl_math.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_errno.h>

#include "third-party/fast-cpp-csv-parser/csv.h"

#define HOPELESS_MANEUVER_ENERGY_COST 100          // Values to return when an interpolation is not available becaue the
#define HOPELESS_MANEUVER_PURSUIT_DURATION 100     // fish probably can't maneuver under the requested conditions.
#define PRINT_AVAILABILITY_WARNINGS false

namespace fs = std::experimental::filesystem;

class ManeuverInterpolation {

private:

    gsl_spline2d *energy_cost_spline, *pursuit_duration_spline;
    gsl_interp_accel *xacc_ec, *yacc_ec, *xacc_pd, *yacc_pd;
    bool ec_spline_initiated = false;
    bool pd_spline_initiated = false;

    bool tables_available = true;   // Turns to false if no precalculated tables are available for this fish size

    double velocity_cms, temperature_C, fork_length_cm;
    fs::path *base_path;         // Folder containing the entire hierarchy of maneuver interpolation tables
    fs::path files_folder();    // Returns the folder containing the specific csvs for this interpolation object
    void load_csv(bool is_energy_cost);

public:

    ManeuverInterpolation(double velocity_cms, double temperature_C, double fork_length_cm, fs::path *base_path);

    virtual ~ManeuverInterpolation();

    double interpolate(double xmr, double ymr, bool is_energy_cost);

};

#endif //DRIFTMODELC_MANEUVERINTERPOLATION_H
