#include <string>
#include <iostream>
#include <cmath>

#include "gradalg.h"
#include "tool_laplace.h"
#include "random_gen.hpp"
#include "inverse_laplace_operator.hpp"

int main(){

  //double * M, * f, * b;
  double * eigenvs;
  double sigma = 0.6,  s = -0.5;
  double r_hat = 1.e-16;

  int i, num_iter, L = 6;

  double * M = new double[L * L];
  double * f = new double[L];
  double * b = new double[L];
  double * check_sol = new double[L];


  std::ofstream standard_implemetation;
  std::ofstream efficient_implemetation;

  /* Solution using Conjugate Gradient */

  init_laplace_matrix(M, sigma, s, L);

  fill_source(b, 2.2, 0.5, L);

  conj_gradient_algorithm(M, f, b, r_hat, L, &num_iter);

  inverse_laplace_operator(check_sol, b, sigma, L, L);


  std::cout <<"\n\t Solution: \t\t Check: " << "\t" << std::endl;

  for(i = 0; i < L; i++){
    std::cout << "\t" << f[i] << "\t\t"
    << check_sol[i] << "\t"
    << std::endl;
  }

  delete [] M;
  delete [] f;
  delete [] b;
  delete [] check_sol;



}
