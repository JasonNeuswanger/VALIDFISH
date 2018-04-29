//
// Created by Jason Neuswanger on 2/1/18.
//

#include "Forager.h"


struct volume_integrand_inner_params {gsl_function *func; double theta;};
double volume_integrand_inner(double rho, void *params) {
    struct volume_integrand_inner_params *p = (struct volume_integrand_inner_params *) params;
    const double dV = gsl_pow_2(rho) * sin(p->theta);
    auto *f3d = (gsl_function_pp_3d<double(double, double *, double *)> *) (p->func);
    return f3d->function(rho, f3d->params) * dV;
}

struct volume_integrand_middle_params {gsl_function *func; double phi; double zmin; double zmax; double radius; double min_rho; double max_rho; };
double volume_integrand_middle(double theta, void *params) {
    struct volume_integrand_middle_params *p = (struct volume_integrand_middle_params *) params;
    struct volume_integrand_inner_params inner_params = {p->func, theta};
    gsl_function F = {&volume_integrand_inner, &inner_params};
    double result, error;
    const double rhoMin = isnan(p->min_rho) ? 0.0 : fmax(0.0, p->min_rho);
    const double zLim = (p->phi > 0.0) ? p->zmax : p->zmin;
    const double rhoLim = (theta != 0.0 && p->phi != 0.0) ? fabs(zLim / (sin(p->phi) * sin(theta))) : p->radius;
    const double rhoMax = isnan(p->max_rho) ? fmin(rhoLim, p->radius) : fmin(rhoLim, p->max_rho);
    auto *f3d = (gsl_function_pp_3d<double(double, double *, double *)> *) (p->func);
    *(f3d->y) = theta; // Using 'y' in gsl_function_pp_3d as a placeholder to pass theta to the inner integral
#if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
        gsl_integration_qags(&F, rhoMin, rhoMax, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
        gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(&F, rhoMin, rhoMax, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    return result;
}

struct volume_integrand_outer_params {gsl_function *func; double zmin; double zmax; double radius; double theta;
                                      double min_rho; double max_rho; double min_theta; double max_theta; };
double volume_integrand_outer(double phi, void *params) {
    struct volume_integrand_outer_params *p = (struct volume_integrand_outer_params *) params;
    struct volume_integrand_middle_params middle_params = {p->func, phi, p->zmin, p->zmax, p->radius, p->min_rho, p->max_rho};
    gsl_function F = {&volume_integrand_middle, &middle_params};
    double result, error;
    const double thetaMin = isnan(p->min_theta) ? 0.0 : p->min_theta / 2.0;
    const double thetaMax = isnan(p->max_theta) ? p->theta / 2.0 : p->max_theta / 2.0;
    auto *f3d = (gsl_function_pp_3d<double(double, double *, double *)> *) (p->func);
    *(f3d->x) = phi; // Using 'x' in gsl_function_pp_3d as a placeholder to pass phi to the inner integrals
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
            gsl_integration_qags(&F, thetaMin, thetaMax, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
            gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(&F, thetaMin, thetaMax, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    return result;
}

double Forager::integrate_over_volume(gsl_function *func, double min_rho, double max_rho, double min_theta, double max_theta) {
    /* Integrates a 3-D function f(rho, theta, phi) over the search volume. Pass NAN for the rho/theta limits for the default
     * behavior of integrating over the full search volume in one or both dimensions. */
    struct volume_integrand_outer_params outer_params = {func, bottom_z, surface_z, max_radius, theta, min_rho, max_rho, min_theta, max_theta};
    gsl_function F = {&volume_integrand_outer, &outer_params};
    double result, error;
    const double phiMin = -M_PI_2;
    const double phiMax = M_PI_2;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
            gsl_integration_qags(&F, phiMin, phiMax, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
            gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(&F, phiMin, phiMax, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    return 2 * result;
}

double Forager::volume_within_radius(double r) {
    double inner_theta, inner_phi;  // placeholders for the inner integrands
    auto one_as_func_3d = [](double rho, double *theta, double *phi)->double{ return 1.0; };
    gsl_function_pp_3d<decltype(one_as_func_3d)> Fintegrand(one_as_func_3d, &inner_theta, &inner_phi);
    gsl_function *F = static_cast<gsl_function*>(&Fintegrand);
    return integrate_over_volume(F, 0, r, NAN, NAN);
}

double Forager::cross_sectional_area() {
    auto area_element = [](double x)->double{ return 1; }; //    auto area_element = [] (double z, void *params) { return z; };
    gsl_function_pp<decltype(area_element)> Fp(area_element);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    return integrate_over_xz_plane(F, true);
}

struct integrate_over_xz_plane_inner_params {gsl_function *F; double rho; double bottom_z; double surface_z; bool integrand_is_1d; };

double integrate_over_xz_plane_inner(double x, void *params) {
    struct integrate_over_xz_plane_inner_params *p = (struct integrate_over_xz_plane_inner_params *) params;
    double result, error;
    const double z_circle_edge = sqrt(pow(p->rho, 2) - pow(x, 2));
    const double zMin = fmax(-z_circle_edge, p->bottom_z);
    const double zMax = fmin(z_circle_edge, p->surface_z);
    if (!p->integrand_is_1d) {
        auto *f3d = (gsl_function_pp_3d<double(double, double *, double *)> *) (p->F);  // as of here, p->F->x has no value... receives it below
        *(f3d->x) = x; // For 3D versions later, if any, change this to *(f3d->y) = y in the y loop
    }
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
        gsl_integration_qags(p->F, zMin, zMax, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
        gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(p->F, zMin, zMax, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    return result;
}

double Forager::integrate_over_xz_plane(gsl_function *func, bool integrand_is_1d) {
    const double rho = (theta >= M_PI) ? max_radius : max_radius * sin(theta / 2.0);
    struct integrate_over_xz_plane_inner_params inner_params = {func, rho, bottom_z, surface_z, integrand_is_1d};
    gsl_function F;
    F.function = &integrate_over_xz_plane_inner;
    F.params = &inner_params;
    double result, error;
    const double xMin = 0.0; // start at 0 because of symmetry, multiply answer by 2 below
    const double xMax = rho;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
        gsl_integration_qags(&F, xMin, xMax, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
        gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(&F, xMin, xMax, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    return 2.0 * result;
}

double Forager::integrate_energy_cost_over_prey_path(double x, double z, std::shared_ptr<PreyType> pt, bool is_energy_cost) {
    const double prey_radius = pt->get_max_attended_distance();
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(prey_radius);
    const double rhosq = (theta < M_PI) ? gsl_pow_2(prey_radius * sin(theta / 2)) : rsq; // rho = max radius in lateral direction
    if (xsq + zsq >= rhosq) { return 0; }    // if (x,z) are outside search volume return 0
    const double y0 = sqrt(rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    const double v = water_velocity(z); // Duplicating the passage_time function within this one because we need v below
    const double T = (y0 - yT) / v;     // and this way there's no need to calculate water_velocity(z) twice.
    auto integrand = [this, x, z, pt, y0, v, prey_radius, is_energy_cost](double t)->double{
        const double tauval = tau(t, x, z, pt);
        const double y = y0 - v * t;
        assert(prey_radius*prey_radius >= x*x + z*z + y*y);
        if (isfinite(tauval)) {
            const double v = (water_velocity(z) + focal_velocity) / 2;  // Calculate cost from avg of focal & prey position velocities
            return maneuver_cost(x, y, z, v, is_energy_cost) * exp(-t / tauval) / tauval;
        } else {
            return 0.;
        }
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
        gsl_integration_qags(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
        gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    if (DIAG_NANCHECKS) assert(isfinite(result));
    return result;
}

double Forager::passage_time(double x, double z, std::shared_ptr<PreyType> pt) {
    /* Returns the duration of time (s) for which the item at (x, z) is passing through the search volume
     * This is only used in integrate_detection_pdf below, but it's also shared in the Python module for analysis. */
    const double prey_radius = pt->get_max_attended_distance();
    const double xsq = gsl_pow_2(x);
    const double zsq = gsl_pow_2(z);
    const double rsq = gsl_pow_2(prey_radius);
    const double rhosq = (theta < M_PI) ? gsl_pow_2(prey_radius * sin(theta / 2)) : rsq; // rho = max prey_radius in lateral direction
    if (xsq + zsq >= rhosq) { return 0; }    // if (x,z) are outside search volume return 0
    const double y0 = sqrt(rsq - xsq - zsq);
    const double yT = fmax(-y0, cot(theta/2) * sqrt(xsq + zsq));
    return (y0 - yT) / water_velocity(z);
}

double Forager::integrate_detection_pdf(double x, double z, std::shared_ptr<PreyType> pt) {
    /* Integrates the detection pdf over the prey's path from t=0 to t=T */
    const double T = passage_time(x, z, pt);
    if (T <= 0) { return 0; }  // If (x, z) are outside the search volume, detection probability is 0.
    auto integrand = [this, x, z, pt](double t)->double{
        double tauval = tau(t, x, z, pt);
        if (isfinite(tauval)) {
            return exp(-t / tauval) / tauval;
        } else {
            return 0.;
        }
    };
    gsl_function_pp<decltype(integrand)> Fp(integrand);
    gsl_function *F = static_cast<gsl_function*>(&Fp);
    double result, error;
    #if USE_ADAPTIVE_INTEGRATION
        gsl_integration_workspace *w = gsl_integration_workspace_alloc(QUAD_SUBINT_LIM);
        gsl_integration_qags(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, QUAD_SUBINT_LIM, w, &result, &error);
        gsl_integration_workspace_free(w);
    #else
        size_t neval;
        gsl_integration_qng(F, 0, T, QUAD_EPSABS, QUAD_EPSREL, &result, &error, &neval);
    #endif
    assert(isfinite(result));
    if (!isfinite(result)) {
        printf("TAU_ERROR Detection PDF integral returned non-finite result %.3f for T=%.3f at (x,z)=(%.3f,%.3f) for prey type %s.\n", result, T, x, z, pt->name.c_str());
    }
    return result;
}

double* Forager::random_xz() {
    /* Returns a randomly selected (x, z) coordinate within the foraging volume of the fish.
     * It's not quite uniformly distributed, but good enough for randomized testing. */
    static double ret[2];
    double rho = (theta < M_PI) ? max_radius * sin(theta / 2) : max_radius;
    double x = rho * (double)rand() / RAND_MAX;
    double zbl = fmax(-sqrt(gsl_pow_2(rho) - gsl_pow_2(x)), bottom_z);
    double zbu = fmin(sqrt(gsl_pow_2(rho) - gsl_pow_2(x)), surface_z);
    double z = zbl + ((double)rand() / RAND_MAX) * (zbu - zbl);
    ret[0] = x;
    ret[1] = z;
    return ret;
}