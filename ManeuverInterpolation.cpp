//
// Created by Jason Neuswanger on 2/19/18.
//

#include "ManeuverInterpolation.h"

ManeuverInterpolation::ManeuverInterpolation(double velocity_cms, double temperature_C, double fork_length_cm, fs::path *base_path) {
    assert(fs::exists(*base_path));
    assert(fs::is_directory(*base_path));
    this->base_path = base_path;
    this->velocity_cms = velocity_cms;
    this->temperature_C = temperature_C;
    this->fork_length_cm = fork_length_cm;
    xacc_ec = gsl_interp_accel_alloc(); // Allocating memory for caches to speed interpolation table lookups
    yacc_ec = gsl_interp_accel_alloc();
    xacc_pd = gsl_interp_accel_alloc();
    yacc_pd = gsl_interp_accel_alloc();
    load_csv(true);
    load_csv(false);
}

ManeuverInterpolation::~ManeuverInterpolation() {
    if (ec_spline_initiated) {
        gsl_spline2d_free(energy_cost_spline);
        gsl_interp_accel_free(xacc_ec);
        gsl_interp_accel_free(yacc_ec);
    }
    if (pd_spline_initiated) {
        gsl_spline2d_free(pursuit_duration_spline);
        gsl_interp_accel_free(xacc_pd);
        gsl_interp_accel_free(yacc_pd);
    }
}

fs::path ManeuverInterpolation::files_folder() {
    // First, find the folder name most closely matching the fish's fork length in cm
    double nearest_fork_length = NAN;
    double nearest_fork_length_distance = INFINITY;
    std::string nearest_fork_length_str;
    for (const auto& entry : fs::directory_iterator(*base_path)) {
        auto filename = entry.path().filename().string();
        if (fs::is_directory(entry.status())) {
            std::string fork_length_str = filename.substr(filename.find("_")+1, filename.length());
            const double folder_fork_length = std::stof(fork_length_str);
            const double fork_length_distance = fabs(fork_length_cm - folder_fork_length);
            if (fork_length_distance < nearest_fork_length_distance) {
                nearest_fork_length_distance = fork_length_distance;
                nearest_fork_length = folder_fork_length;
                nearest_fork_length_str = fork_length_str;
            }
        }
    }
    assert(!isnan(nearest_fork_length));
    fs::path fork_length_path = *base_path / ("fl_" + nearest_fork_length_str);
    assert(fs::exists(fork_length_path));
    assert(fs::is_directory(fork_length_path));
    if (nearest_fork_length_distance > 1.1) {
        tables_available = false;
        if (PRINT_AVAILABILITY_WARNINGS) {
            fprintf(stderr, "WARNING: Maneuver tables not available for fork length %.1f cm. ", fork_length_cm);
            fprintf(stderr, "Nearest is %.1f cm.\n", nearest_fork_length);
        }
    }
    // Second, find the folder name most closely matching the intended velocity for the maneuver in cm/s
    double nearest_velocity = NAN;
    double nearest_velocity_distance = INFINITY;
    std::string nearest_velocity_str;
    for (const auto& entry : fs::directory_iterator(fork_length_path)) {
        auto filename = entry.path().filename().string();
        if (fs::is_directory(entry.status())) {
            std::string velocity_str = filename.substr(filename.find("_")+1, filename.length());
            const double folder_velocity = std::stof(velocity_str);
            const double velocity_distance = fabs(velocity_cms - folder_velocity);
            if (velocity_distance < nearest_velocity_distance) {
                nearest_velocity_distance = velocity_distance;
                nearest_velocity = folder_velocity;
                nearest_velocity_str = velocity_str;
            }
        }
    }
    assert(!isnan(nearest_velocity));
    fs::path velocity_path = fork_length_path / ("fcs_" + nearest_velocity_str);
    assert(fs::exists(velocity_path));
    assert(fs::is_directory(velocity_path));
    if (nearest_velocity_distance > 2.6) {
        tables_available = false;
        if (PRINT_AVAILABILITY_WARNINGS) {
            fprintf(stderr, "WARNING: Maneuver tables not available for fork length %.1f cm ", nearest_fork_length);
            fprintf(stderr, "at velocity %.1f cm/s. ", velocity_cms);
            fprintf(stderr, "Nearest is %.1f cm/s instead.\n", nearest_velocity);
        }
    }
    // Third, find the folder name most closely matching the fish's temperature
    double nearest_temperature = NAN;
    double nearest_temperature_distance = INFINITY;
    std::string nearest_temperature_str;
    for (const auto& entry : fs::directory_iterator(velocity_path)) {
        auto filename = entry.path().filename().string();
        if (fs::is_directory(entry.status())) {
            std::string temperature_str = filename.substr(filename.find("_")+1, filename.length());
            const double folder_temperature = std::stof(temperature_str);
            const double temperature_distance = fabs(temperature_C - folder_temperature);
            if (temperature_distance < nearest_temperature_distance) {
                nearest_temperature_distance = temperature_distance;
                nearest_temperature = folder_temperature;
                nearest_temperature_str = temperature_str;
            }
        }
    }
    assert(!isnan(nearest_temperature));
    fs::path temperature_path = velocity_path / ("t_" + nearest_temperature_str);
    assert(fs::exists(temperature_path));
    assert(fs::is_directory(temperature_path));
    if (nearest_temperature_distance > 1.1) {
        tables_available = false;
        if (PRINT_AVAILABILITY_WARNINGS) {
            fprintf(stderr, "WARNING: Maneuver tables not available for fork length %.1f cm ", nearest_fork_length);
            fprintf(stderr, "at velocity %.1f cm/s and temperature %.1f C. ", nearest_velocity, temperature_C);
            fprintf(stderr, "Nearest is %.1f C.\n", nearest_temperature);
        }
    }
    return temperature_path.string();
}

std::string ManeuverInterpolation::source_filename(bool is_energy_cost) {
    fs::path file_name = (is_energy_cost) ? fs::path("energy_cost.csv") : fs::path("pursuit_duration.csv");
    fs::path file_path = files_folder() / file_name;
    return file_path.string();
}

void ManeuverInterpolation::load_csv(bool is_energy_cost) {
    // Create the LineReader object
    fs::path file_name = (is_energy_cost ? fs::path("energy_cost.csv") : fs::path("pursuit_duration.csv"));
    fs::path file_path = files_folder() / file_name;
    if (!tables_available) { return; } // Don't even bother building structures if tables aren't available for the conditions
    io::LineReader in(file_path);
    // Read the first row (i.e., first line) and count the number of commas to figure out how many columns there are
    char *xline = in.next_line();
    std::string xline_str(xline);
    const size_t ncols = (size_t) std::count(xline_str.begin(), xline_str.end(), ',') - 1;
    // Process the contents of the first row into the x vector
    std::vector<double> x, y;
    io::detail::custom_parse_line_doubles< io::trim_chars<' ', '\t'> , io::no_quote_escape<','>>(xline, x, ncols+1);
    x.erase(x.begin()); // remove the first (nan) element from the x values
    // Read the remaining rows into a vector of vectors
    std::vector<std::vector<double>> rows;
    while(char *line = in.next_line()){
        std::vector<double> values;
        io::detail::custom_parse_line_doubles< io::trim_chars<' ', '\t'> , io::no_quote_escape<','>>(line, values, ncols+1);
        rows.emplace_back(values);
    }
    const size_t nrows = rows.size(); // number of data rows (header excluded by being read separately)
    // Go through each row and pull out the y value into a separate vector
    for (auto&row : rows) {
        y.emplace_back(row.at(0));
        row.erase(row.begin());
    }
    // Copy the x and y vectors into c arrays usable by gsl by copying from the arrays they use internally
    double* xraw = &x[0];
    double* yraw = &y[0];
    double xa[ncols];
    double ya[nrows];
    memcpy(xa, xraw, ncols * sizeof(xraw));
    memcpy(ya, yraw, nrows * sizeof(yraw));
    // Initialize the gsl spline-related objects
    auto *za = (double *) malloc(ncols * nrows * sizeof(double));
    const gsl_interp2d_type *T = gsl_interp2d_bicubic;
    gsl_spline2d **spline;  // pointer to one of the instance variable spline objects
    bool *spline_initiated;
    if (is_energy_cost) {
        spline = &energy_cost_spline;
        spline_initiated = &ec_spline_initiated;
    } else {
        spline = &pursuit_duration_spline;
        spline_initiated = &pd_spline_initiated;
    }
    if (!*spline_initiated) {
        *spline = gsl_spline2d_alloc(T, ncols, nrows);
        *spline_initiated = true;
    }
    for (size_t i = 0; i < ncols; ++i) {
        for (size_t j = 0; j < nrows; ++j) {
            gsl_spline2d_set(*spline, za, i, j, rows.at(j).at(i));
        }
    }
    gsl_spline2d_init(*spline, xa, ya, za, ncols, nrows);
    free(za);
}

double ManeuverInterpolation::interpolate(double xmr, double ymr, bool is_energy_cost) {
    /* Takes inputs in maneuver model coordinates, i.e. in the +y half of x-y plane where +x is the downstream
     * direction and +y is lateral, with the positive direction on the right when facing upstream. Units for the
     * spatial inputs are centimeters. */
    if (!tables_available) {
        return (is_energy_cost) ? HOPELESS_MANEUVER_ENERGY_COST : HOPELESS_MANEUVER_PURSUIT_DURATION;
    }
    gsl_spline2d **spline;
    gsl_interp_accel **xacc, **yacc;
    if (is_energy_cost) {
        spline = &energy_cost_spline;
        xacc = &xacc_ec;
        yacc = &yacc_ec;
    } else {
        spline = &pursuit_duration_spline;
        xacc = &xacc_pd;
        yacc = &yacc_pd;
    }
    double mc;
    int status = gsl_spline2d_eval_e(*spline, ymr, xmr, *xacc, *yacc, &mc);
    if (status == GSL_SUCCESS) {
        if (isfinite(mc)) {
            return mc;
        } else {
            if (PRINT_AVAILABILITY_WARNINGS) {
                fprintf(stderr, "WARNING: Requested a maneuver cost outside the interpolated domain. In maneuver plane coordinates, (xmr, ymr) = (%.4f, %.4f) cm. Returning 10 s or J/s.\n", xmr, ymr);
            }
            return (is_energy_cost) ? HOPELESS_MANEUVER_ENERGY_COST : HOPELESS_MANEUVER_PURSUIT_DURATION;
        }
    } else if (status == GSL_EDOM) {
        if (PRINT_AVAILABILITY_WARNINGS) {
            fprintf(stderr, "WARNING: Requested a maneuver cost outside the interpolated domain. In maneuver plane coordinates, (xmr, ymr) = (%.4f, %.4f) cm. Returning 10 s or J/s.\n", xmr, ymr);
        }
        return (is_energy_cost) ? HOPELESS_MANEUVER_ENERGY_COST : HOPELESS_MANEUVER_PURSUIT_DURATION;
    } else {
        fprintf(stderr, "ERROR: GSL error code %d in maneuver cost interpolation. Returning 10 s or 10 J/s instead.\n", status);
        return (is_energy_cost) ? HOPELESS_MANEUVER_ENERGY_COST : HOPELESS_MANEUVER_PURSUIT_DURATION;
    }
}