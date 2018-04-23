//
// Created by Jason Neuswanger on 2/2/18.
//

#ifndef DRIFTMODELC_UTILITY_H
#define DRIFTMODELC_UTILITY_H

#include <bitset>
#include <cmath>
#include "PreyType.h"
#include "Forager.h"
#include <gsl/gsl_errno.h>
#include <random>
#include <numeric>



double cot(double x);

long long xzpciec_hash_key(double x, double z, PreyType *pc, bool is_energy_cost);

void print_gsl_errors(const char * reason, const char * file, int line, int gsl_errno);
void ignore_gsl_errors(const char * reason, const char * file, int line, int gsl_errno);

struct cartesian_3D_coords { double x; double y; double z; };

cartesian_3D_coords cartesian_from_spherical(double rho, double theta, double phi);

double fRand(double fMin, double fMax);

double trim_to_bounds(double value, std::array<double,2>bounds);

#endif //DRIFTMODELC_UTILITY_H
