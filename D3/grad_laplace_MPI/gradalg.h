#ifndef __GRAD__
#define __GRAD__

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>



void gradient_algorithm(double * A, double * x, double * b, double r_hat_targ,
  int N, int * num_iter);

void conj_gradient_algorithm(double * A, double * x, double * b, double r_hat_targ,
  int N, int * num_iter);

void errors_conj_grad_alg(double * A, double * x, double * b, double r_hat_targ,
  int N);

void efficient_conj_grad_alg(double * x, double * b,
    double sigma, double s, double r_hat_targ, int N, int * num_iter);

void efficient_conj_grad_alg_MPI(double * x, double * b,
    double sigma, double s, double r_hat_targ, int N, int * num_iter, int rank, int size);


#endif
