//
// Created by Jason Neuswanger on 2/16/18.
//

#include "Swimmer.h"

static const double smr_params[26][5] = {{45.5499,-1.6385,-0.3505,0.0251,-0.0005},{49.6992,-3.7997,0.2781,-0.0489,0.0025},{53.8485,-5.9608,0.9068,-0.1228,0.0055},{57.9978,-8.1219,1.5355,-0.1968,0.0084},{62.1471,-10.283,2.1641,-0.2708,0.0114},{66.2965,-12.4441,2.7928,-0.3448,0.0144},{71.9988,-12.783,2.6376,-0.321,0.0135},{77.7011,-13.1218,2.4824,-0.2972,0.0126},{83.4034,-13.4606,2.3272,-0.2735,0.0117},{89.1057,-13.7995,2.172,-0.2497,0.0108},{94.808,-14.1383,2.0168,-0.2259,0.0099},{103.716,-13.7171,1.3149,-0.1193,0.0051},{112.623,-13.2959,0.613,-0.0127,0.0002},{121.531,-12.8747,-0.0889,0.0938,-0.0046},{130.439,-12.4535,-0.7908,0.2004,-0.0094},{139.346,-12.0323,-1.4927,0.307,-0.0142},{151.881,-12.9877,-1.682,0.3381,-0.0156},{164.415,-13.9431,-1.8713,0.3692,-0.017},{176.95,-14.8985,-2.0607,0.4003,-0.0184},{189.484,-15.8539,-2.25,0.4314,-0.0198},{202.019,-16.8093,-2.4393,0.4625,-0.0212},{218.934,-19.3705,-1.6324,0.2911,-0.0107},{235.848,-21.9317,-0.8255,0.1196,-0.0002},{252.763,-24.4929,-0.0186,-0.0518,0.0103},{269.678,-27.0541,0.7884,-0.2233,0.0207},{286.592,-29.6153,1.5953,-0.3947,0.0312}};

static const double amr_params[26][5] = {{439.744,0.835,-0.5782,0.1335,-0.0096},{450.73,1.1548,-0.9083,0.1983,-0.0136},{461.716,1.4747,-1.2384,0.263,-0.0176},{472.702,1.7945,-1.5684,0.3278,-0.0216},{483.688,2.1143,-1.8985,0.3925,-0.0255},{494.674,2.4341,-2.2286,0.4573,-0.0295},{531.434,2.555,-1.8418,0.2961,-0.0199},{568.194,2.6759,-1.455,0.1349,-0.0102},{604.953,2.7968,-1.0683,-0.0263,-0.0006},{641.713,2.9177,-0.6815,-0.1874,0.009},{678.473,3.0387,-0.2948,-0.3486,0.0187},{741.436,7.4827,-6.8461,0.7548,-0.0356},{804.399,11.9268,-13.3975,1.8581,-0.0899},{867.362,16.3708,-19.9489,2.9615,-0.1442},{930.325,20.8149,-26.5003,4.0649,-0.1985},{993.288,25.2589,-33.0516,5.1682,-0.2528},{977.569,20.0484,-27.7516,4.2329,-0.2033},{961.851,14.8378,-22.4515,3.2974,-0.1538},{946.133,9.6273,-17.1515,2.3621,-0.1044},{930.415,4.4168,-11.8514,1.4267,-0.0549},{914.697,-0.7937,-6.5514,0.4913,-0.0054},{906.587,-0.927,-6.8678,0.6507,-0.0175},{898.478,-1.0602,-7.1842,0.8101,-0.0295},{890.368,-1.1935,-7.5006,0.9696,-0.0416},{882.259,-1.3268,-7.8171,1.129,-0.0536},{874.149,-1.46,-8.1335,1.2884,-0.0657}};

static const double u_ms_params[26][5] = {{16.0807,9.7663,-0.9914,0.1493,0.0005},{16.8226,11.0578,-1.8682,0.3386,-0.0115},{17.5645,12.3492,-2.7451,0.5279,-0.0235},{18.3064,13.6407,-3.622,0.7173,-0.0355},{19.0483,14.9322,-4.4989,0.9066,-0.0475},{19.7902,16.2237,-5.3757,1.0959,-0.0596},{21.7221,14.583,-4.4013,0.9199,-0.0487},{23.6541,12.9423,-3.427,0.7439,-0.0378},{25.5861,11.3016,-2.4526,0.568,-0.027},{27.5181,9.6609,-1.4782,0.392,-0.0161},{29.4501,8.0202,-0.5038,0.216,-0.0053},{31.2651,7.6782,-0.109,0.1486,-0.0008},{33.0801,7.3363,0.2857,0.0813,0.0036},{34.8951,6.9943,0.6804,0.014,0.008},{36.7101,6.6524,1.0751,-0.0534,0.0125},{38.5251,6.3104,1.4699,-0.1207,0.0169},{38.0136,6.759,1.1385,-0.0551,0.0127},{37.5021,7.2076,0.8072,0.0106,0.0085},{36.9906,7.6561,0.4759,0.0763,0.0043},{36.4791,8.1047,0.1445,0.1419,0.0001},{35.9676,8.5533,-0.1868,0.2076,-0.004},{35.2118,9.1098,-0.4654,0.251,-0.0064},{34.4561,9.6662,-0.7441,0.2945,-0.0088},{33.7004,10.2227,-1.0227,0.3379,-0.0112},{32.9446,10.7792,-1.3014,0.3814,-0.0135},{32.1889,11.3357,-1.58,0.4248,-0.0159}};

Swimmer::Swimmer(double fork_length_cm,
                 double mass_g,
                 unsigned temperature_C,
                 std::string *maneuver_interpolation_csv_base_path) {
    if (fork_length_cm == -1) { // Pass a fork length of -1 to infer it instead
        double total_length_mm = pow(10, (0.3307 * (5.023 + log10(mass_g))));
        double fork_length_mm = (total_length_mm + 0.027) / 1.072;
        this->fork_length_cm = fork_length_mm / 10;
    } else {
        this->fork_length_cm = fork_length_cm;
    }
    this->mass_g = mass_g;
    this->temperature_C = temperature_C;
    this->maneuver_interpolation_csv_base_path = fs::path(maneuver_interpolation_csv_base_path->c_str());
    process_parameter_updates(true);
}

Swimmer::Swimmer(Swimmer *otherSwimmer) {
    // Deep copy an existing Swimmer, minus maneuver splines/accelerators -- they'll be rebuilt as needed.
    mass_g = otherSwimmer->mass_g;
    fork_length_cm = otherSwimmer->fork_length_cm;
    temperature_C = otherSwimmer->temperature_C;
    maneuver_interpolation_csv_base_path = otherSwimmer->maneuver_interpolation_csv_base_path;
    process_parameter_updates(true);
}

void Swimmer::modify_parameters(double fork_length_cm,
                                double mass_g,
                                unsigned temperature_C) {
    assert(fork_length_cm > 0);
    assert(mass_g > 0);
    assert(temperature_C >= 0);
    bool reset_interpolations = ((fork_length_cm != this->fork_length_cm) || (temperature_C != this->temperature_C));
    this->fork_length_cm = fork_length_cm;
    this->mass_g = mass_g;
    this->temperature_C = temperature_C;
    process_parameter_updates(reset_interpolations);
}

void Swimmer::process_parameter_updates(bool reset_interpolations) {
    SMR = brett_glass_regression_value(0);   // standard metabolic rate (mgO2 * kg/h)
    AMR = brett_glass_regression_value(1);   // active metabolic rate (mgO2 * kg/h), i.e. oxygen consumption at u_ms
    u_ms = 0.01 * brett_glass_regression_value(2);  // maximum sustainable swimming speed in m/s, converted from cm/s
    if (reset_interpolations) { interps.erase(interps.begin(), interps.end()); }
}

double Swimmer::brett_glass_regression_value(unsigned which_params) {
    const double (*params)[26][5];
    switch (which_params) {
        case 0:
            params = &smr_params;
            break;
        case 1:
            params = &amr_params;
            break;
        case 2:
            params = &u_ms_params;
            break;
        default:
            throw std::runtime_error("Invalid brett_glass params specification.");
    }
    const double *b = (*params)[temperature_C];
    const double b1 = b[0];
    const double b2 = b[1];
    const double b3 = b[2];
    const double b4 = b[3];
    const double b5 = b[4];
    double logmass = log((mass_g < 0.8) ? 0.8 : mass_g); // round mass up to 0.8 for tiny fish to avoid regression problems
    return b1 + b2 * logmass + b3 * gsl_pow_2(logmass) + b4 * gsl_pow_3(logmass) + b5 * gsl_pow_4(logmass);
}

double Swimmer::steady_swimming_cost(double swimming_speed) {
    return (1/3600.0) * (mass_g / 1000) * 14.06 * SMR * pow((AMR/SMR),(swimming_speed/u_ms));
}

double Swimmer::maneuver_cost(double x, double y, double z, double v, bool is_energy_cost) {
    const double xmr = -100 * y;                                  // transforming input coordinates into maneuver
    const double ymr = 100 * sqrt(gsl_pow_2(x) + gsl_pow_2(z));   // model coordinates and units from m to cm
    double interp_velocity = round((v * 100) - (((int) round(v * 100) - 1) % 5)); // Interpolate at velocity matching 1, 6, 11, 16, 21, 26, 31, 36 .... cm/s
    if (interps.count(interp_velocity) == 0) {  // If interpolations for this velocity aren't loaded yet, load them
        interps.insert(std::pair<double, std::shared_ptr<ManeuverInterpolation>>(interp_velocity, new ManeuverInterpolation(interp_velocity, temperature_C, fork_length_cm, &maneuver_interpolation_csv_base_path)));
    }
    double result = interps.at(interp_velocity)->interpolate(xmr, ymr, is_energy_cost);
    assert(isfinite(result));
//    if (result <= 0) {
//        std::string cost_type_label = (is_energy_cost) ? "energy" : "pursuit duration";
//        std::string file_path = interps.at(interp_velocity)->source_filename(is_energy_cost);
//        printf("Interpolated a nonpositive maneuver cost of %.6f J at (x, y, z) = (%.4f, %.4f, %.4f) with v=%.6f and cost type %s.\n", result, x, y, z, v, cost_type_label.c_str());
//        printf("In maneuver model coordinates, (xmr, ymr) = (%.4f, %.4f) and interpolation file is %s.\n", xmr, ymr, file_path.c_str());
//    }
//    assert(result > 0);
    return fmax(result, 1e-6);  // todo Figure out why I was sometimes getting negative costs interpolated from very small positive costs and fix it. Plot stuff in Python.
}