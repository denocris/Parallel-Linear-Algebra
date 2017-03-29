#ifndef __TOOLS__
#define __TOOLS__

#include <string>
#include <iostream>
#include <cmath>

double scalar_prod(double *, double *, int);

double * mat_vec_prod(const double * A, const double * x, const int N);

double vector_norm(const double * x, const int N);

#endif
