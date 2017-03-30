#ifndef __TOOLS__
#define __TOOLS__

#include <string>
#include <iostream>
#include <cmath>

#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

void init_laplace_matrix(double* M, double sigma, double s, int L);

double seconds();

void compute_eigenvalues(double* eigenv, double sigma, int L);

#endif
