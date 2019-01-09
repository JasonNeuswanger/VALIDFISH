//
// Created by Jason Neuswanger on 5/15/18.
//

#include "Forager.h"

void Forager::set_strategy_bounds() {
    strategy_bounds[s_sigma_A][0] = 0.1;
    strategy_bounds_notes[s_sigma_A][0] = "WARNING: Unrealistic forward attention concentration of spatial attention";
    strategy_bounds[s_sigma_A][1] = 4.0;
    strategy_bounds_notes[s_sigma_A][0] = "Approximately uniform spatial attention distribution";

    strategy_bounds[s_mean_column_velocity][0] = 0.01;
    strategy_bounds_notes[s_mean_column_velocity][0] = "WARNING: Selected minimum water velocity, effectively stillwater";
    strategy_bounds[s_mean_column_velocity][1] = 10 * (0.01 * fork_length_cm);
    strategy_bounds_notes[s_mean_column_velocity][1] = "WARNING: Selected arbitrary maximum mean column velocity";

    strategy_bounds[s_inspection_time][0] = 0.167; // physiological max 180 degrees/sec, assume 30 degree average, so 0.167 s just to move the eyes, let alone fixation https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5050213/
    strategy_bounds_notes[s_inspection_time][0] = "Inspection time equal to fastest possible average saccade";
    strategy_bounds[s_inspection_time][1] = 1.0;
    strategy_bounds_notes[s_inspection_time][0] = "WARNING: Selected arbitrary maximum allowed inspection time";

    strategy_bounds[s_discrimination_threshold][0] = -6;
    strategy_bounds_notes[s_discrimination_threshold][0] = "WARNING: Accepts all prey without discrimination";
    strategy_bounds[s_discrimination_threshold][1] = 6; // 6 intrinsic standard deviations beyond the mean preyishness of debris should be plenty
    strategy_bounds_notes[s_discrimination_threshold][1] = "WARNING: Overly discriminating about prey vs debris";

    strategy_bounds[s_search_image][0] = -1;    // Search image value of of -1 to 0 means no search image
    strategy_bounds[s_search_image][1] = 1;     // Search image values between 0 and 1 are split evenly across search-image-eligible categories ordered by number
}

void Forager::set_parameter_bounds() {
    /* All parameter bounds were set via a priori consideration of their intended functions and plausible values in those roles.
     * They therefore cannot be grossly overfitted to serve some hidden meaning very different from the one intended.
     *
     * NOTE ON LOG TRANSFORMS: The parameters are represented in their main variables as they're used in the model. The
     * transform guidance is only for calibration (so the space of possible values can be covered appropriately) and
     * and printing (to see when values are close to the boundary).
     * */

    parameter_bounds[p_delta_0][0] = 1e-6;
    parameter_bounds_notes[p_delta_0][0] = "WARNING: Angular size maximally affects tau";
    parameter_bounds[p_delta_0][1] = 1e-2;
    parameter_bounds_notes[p_delta_0][1] = "WARNING: Angular size negligibly affects tau";
    parameter_log_transforms[p_delta_0] = true;

    parameter_bounds[p_alpha_tau][0] = 1;
    parameter_bounds_notes[p_alpha_tau][0] = "";//"WARNING: Search image negligibly affects tau";
    parameter_bounds[p_alpha_tau][1] = 100;
    parameter_bounds_notes[p_alpha_tau][1] = "WARNING: Search image maximally affects tau";
    parameter_log_transforms[p_alpha_tau] = true;

    parameter_bounds[p_alpha_d][0] = 1;
    parameter_bounds_notes[p_alpha_d][0] = "";// "WARNING: Search image negligibly affects perceptual variance";
    parameter_bounds[p_alpha_d][1] = 100;
    parameter_bounds_notes[p_alpha_d][1] = "WARNING: Search image maximally affects on perceptual variance";
    parameter_log_transforms[p_alpha_d] = true;

    parameter_bounds[p_beta][0] = 0;
    parameter_bounds_notes[p_beta][0] = "WARNING: Set size * inspection time negligibly affects tau";
    parameter_bounds[p_beta][1] = 2;
    parameter_bounds_notes[p_beta][1] = "WARNING: Set size * inspection time maximally affects tau";
    parameter_log_transforms[p_beta] = false;

    parameter_bounds[p_A_0][0] = 0;
    parameter_bounds_notes[p_A_0][0] = "WARNING: Spatial attention negligibly affects tau";
    parameter_bounds[p_A_0][1] = 2;
    parameter_bounds_notes[p_A_0][1] = "WARNING: Spatial attention maximally affects tau";
    parameter_log_transforms[p_A_0] = false;

    parameter_bounds[p_flicker_frequency][0] = 10;
    parameter_bounds_notes[p_flicker_frequency][0] = "WARNING: Flicker frequency set to minimum allowed";
    parameter_bounds[p_flicker_frequency][1] = 70;
    parameter_bounds_notes[p_flicker_frequency][1] = "WARNING: Flicker frequency set to maximum allowed";
    parameter_log_transforms[p_flicker_frequency] = false;

    parameter_bounds[p_tau_0][0] = 1e-2;
    parameter_bounds_notes[p_tau_0][0] = "WARNING: Base tau at minimum allowed value";
    parameter_bounds[p_tau_0][1] = 1e2;
    parameter_bounds_notes[p_tau_0][1] = "WARNING: Base tau at maximum allowed value";
    parameter_log_transforms[p_tau_0] = true;

    parameter_bounds[p_nu_0][0] = 1e-4;
    parameter_bounds_notes[p_nu_0][0] = "WARNING: Loom maximally effects tau";
    parameter_bounds[p_nu_0][1] = 1e-1;
    parameter_bounds_notes[p_nu_0][0] = "WARNING: Loom negligibly effects tau";
    parameter_log_transforms[p_nu_0] = true;

    parameter_bounds[p_discriminability][0] = 1.5;  // Difference in mean preyishness between prey and debris, in units of the (equal) standard deviation of each.
    parameter_bounds_notes[p_discriminability][0] = "WARNING: Prey-debris discrimination is maximally difficult";
    parameter_bounds[p_discriminability][1] = 3.5;
    parameter_bounds_notes[p_discriminability][1] = "WARNING: Prey-debris discrimination is maximally easy";
    parameter_log_transforms[p_discriminability] = false;

    parameter_bounds[p_delta_p][0] = 1e-4;
    parameter_bounds_notes[p_delta_p][0] = "WARNING: Angular size maximally affects perceptual variance";
    parameter_bounds[p_delta_p][1] = 1;
    parameter_bounds_notes[p_delta_p][1] = "WARNING: Angular size negligibly affects perceptual variance";
    parameter_log_transforms[p_delta_p] = true;

    parameter_bounds[p_omega_p][0] = 0;
    parameter_bounds_notes[p_omega_p][0] = "WARNING: Angular velocity negligibly affects perceptual variance";
    parameter_bounds[p_omega_p][1] = 2;
    parameter_bounds_notes[p_omega_p][1] = "WARNING: Angular velocity maximally affects perceptual variance";
    parameter_log_transforms[p_omega_p] = false;

    parameter_bounds[p_ti_p][0] = 1e-2;       //      log scale     -- Scales effect of inspection time on perceptual variance
    parameter_bounds_notes[p_ti_p][0] = "WARNING: Inspection time maximally affects perceptual variance";
    parameter_bounds[p_ti_p][1] = 1e1;
    parameter_bounds_notes[p_ti_p][1] = "WARNING: Inspection time negligibly affects perceptual variance";
    parameter_log_transforms[p_ti_p] = true;

    parameter_bounds[p_sigma_p_0][0] = 1e-2;       //   log scale, really no idea          -- Base constant on which other effects on perceptual variance are multiplied
    parameter_bounds_notes[p_sigma_p_0][0] = "WARNING: Base perceptual variance set to minimum allowed";
    parameter_bounds[p_sigma_p_0][1] = 1e2;
    parameter_bounds_notes[p_sigma_p_0][1] = "WARNING: Base perceptual variance set to maximum allowed";
    parameter_log_transforms[p_sigma_p_0] = true;

}

double Forager::validate(Strategy s, double v) {
    if (v >= strategy_bounds[s][0] && v <= strategy_bounds[s][1]) {
        return v;
    } else if (v < strategy_bounds[s][0]) {
        printf("WARNING! Strategy %s got a value %.8f, which is below the minimum %.8f. Setting to the minimum instead.\n", strategy_names[s].c_str(), v, strategy_bounds[s][0]);
        return strategy_bounds[s][0];
    } else {
        printf("WARNING! Strategy %s got a value %.8f, which is above the maximum %.8f. Setting to the maximum instead.\n", strategy_names[s].c_str(), v, strategy_bounds[s][1]);
        return strategy_bounds[s][1];
    }
}

double Forager::validate(Parameter p, double v) {
    if (v >= parameter_bounds[p][0] && v <= parameter_bounds[p][1]) {
        return v;
    } else if (v < parameter_bounds[p][0]) {
        printf("WARNING! Parameter %s got a value %.8f, which is below the minimum %.8f. Setting to the minimum instead.\n", parameter_names[p].c_str(), v, parameter_bounds[p][0]);
        return parameter_bounds[p][0];
    } else {
        printf("WARNING! Parameter %s got a value %.8f, which is above the maximum %.8f. Setting to the maximum instead.\n", parameter_names[p].c_str(), v, parameter_bounds[p][1]);
        return parameter_bounds[p][1];
    }
}

void Forager::set_single_strategy_bounds(Strategy strategy, double lower_bound, double upper_bound) {
    strategy_bounds[strategy][0] = lower_bound;
    strategy_bounds[strategy][1] = upper_bound;
}

void Forager::fix_single_strategy_bound(Strategy strategy, double fixed_value) {
    strategy_bounds[strategy][0] = fixed_value;
    strategy_bounds[strategy][1] = fixed_value;
}
