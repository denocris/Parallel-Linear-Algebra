#ifndef __GRAD__
#define __GRAD__

#include <string>
#include <iostream>
#include <cmath>



void gradient_algorithm(double * A, double * x, double * b, double r_hat_targ,
  int N, int * num_iter);

void conj_gradient_algorithm(double * A, double * x, double * b, double r_hat_targ,
  int N, int * num_iter);


#endif
