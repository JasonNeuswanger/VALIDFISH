// Wrapper to make it possible to pass lambda functions with captured objects to GSL for integration.
//
// Source of this wrapper: https://stackoverflow.com/questions/13289311/c-function-pointers-with-c11-lambdas/18413206#18413206
// Found via: https://stackoverflow.com/questions/26626158/numerical-integration-of-lambda-function-with-gsl
// Also relevant: https://stackoverflow.com/questions/19168492/using-gsl-functions-defined-in-a-structure

#ifndef DRIFTMODELC_GSL_FUNCTION_PP_H
#define DRIFTMODELC_GSL_FUNCTION_PP_H

#include <gsl/gsl_math.h>

template< typename F >
class gsl_function_pp : public gsl_function {
public:
    explicit gsl_function_pp(const F& func) : _func(func) {
        function = &gsl_function_pp::invoke; // gives (params)->_func(x) but in a different format,
        params=this; // where _func is a pointer to a function of velocity_xz's type (initialized as velocity_xz)
    }
private:
    const F& _func;
    static double invoke(double x, void *params) {
        return static_cast<gsl_function_pp*>(params)->_func(x); // (params) here is "this", so we're saying call this->velocity_xz(x)
    }
};

// Usage of the above for functions that only depend on one variable (including 2-D or 3-D integrations of integrands
// that only depend on one of the variables being integrated over):
//     gsl_function_pp<decltype(velocity_xz)> Fp(velocity_xz);
//     gsl_function *F = static_cast<gsl_function*>(&Fp);


// Version for the inner integrand of a 3D integral, lambdas are f(z, *y, *x) -- z gets passed by the
// inner quadpack integration, x and y are pointers set manually on the way (one or the other may be unused)
// Modified from the above by Jason Neuswanger on 2-4-2018.

template< typename F >
class gsl_function_pp_3d : public gsl_function {
public:
    double *x, *y;
    gsl_function_pp_3d(const F& func, double *y_in, double *x_in) : _func(func) {
        function = &gsl_function_pp_3d::invoke;
        x = x_in;
        y = y_in;
        gfp3dp = {this, y, x};
        params=&gfp3dp;
    }
private:
    struct gsl_function_pp_3d_params {gsl_function_pp_3d *fn; double *y; double *x; } gfp3dp;
    const F& _func;
    static double invoke(double z, void *params) {  // has to accept only x, *params
        auto *p = (struct gsl_function_pp_3d_params *) params;
        //return static_cast<gsl_function_pp_3d*>(p->fn)->_func(z, p->y, p->x); // no longer need the cast, don't know why
        return (p->fn)->_func(z, p->y, p->x);
    }
};

#endif //DRIFTMODELC_GSL_FUNCTION_PP_H
