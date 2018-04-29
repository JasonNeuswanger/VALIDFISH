//
// Created by Jason Neuswanger on 2/2/18.
//

#include "utility.h"

double cot(double x) {
    double s = sin(x);
    double c = cos(x);
    return c / s;
}

long long xzpciec_hash_key(double x, double z, std::shared_ptr<PreyType> pt, bool is_energy_cost) {
    /* This function works by getting 16-bit binary integer representations of each input parameter and
     * concatenating them into a single 64-bit (long long) binary representation of a new integer for the key. */
    long long intx = (long) round(x / MEMOIZATION_PRECISION);
    long long intz = (long) round(z / MEMOIZATION_PRECISION);
    long long ptid = pt->uniqueid;
    long long hash_key_value = (intx << 48) + (intz << 32) + (ptid << 16) + is_energy_cost;
//    std::cout << "val for intx=10k would be " << std::bitset<64>(50) << std::endl;
//    std::cout << "intx = " << intx << " (" << x << ")    intx = " << std::bitset<64>(intx)  << std::endl;
//    std::cout << "intx = " << intz << " (" << z << ")    intx = " << std::bitset<64>(intz)  << std::endl;
//    std::cout << "pcid = " << pcid << "           pcid = " << std::bitset<64>(pcid)  << std::endl;
//    std::cout << "isec = " << is_energy_cost << "           isec = " << std::bitset<64>(is_energy_cost)  << std::endl;
//    std::cout << "hkvl = " << hash_key_value << " = " << std::bitset<64>(hash_key_value)  << std::endl;
    return hash_key_value;
}

void print_gsl_errors(const char * reason, const char * file, int line, int gsl_errno) {
    printf("GSL error code %d (%s) from line %d of file %s.\n", gsl_errno, reason, line, file);
}

void ignore_gsl_errors(const char * reason, const char * file, int line, int gsl_errno) {};

cartesian_3D_coords cartesian_from_spherical(double rho, double theta, double phi) {
    const double sintheta = sin(theta); // because we use it twice
    const double x = rho * sintheta * cos(phi);
    const double y = rho * cos(theta);
    const double z = rho * sintheta * sin(phi);
    return cartesian_3D_coords {x, y, z};
}

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

double trim_to_bounds(double value, std::array<double,2>bounds) {
    if (bounds[0] <= value && value <= bounds[1]) {
        return value;
    } else if (value < bounds[0]) {
        return bounds[0];
    } else {
        return bounds[1];
    }
}