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
  double r_hat_targ = 1.e-16;
  double t_start, t_end;

  int i, num_iter, L = 6;
  int num_rep = 100;
  int L_in = 10, L_fin = 500;

  double * M = new double[L * L];
  double * f = new double[L];
  double * b = new double[L];
  double * check_sol = new double[L];


  std::ofstream standard_implemetation;
  std::ofstream efficient_implemetation;

  /* Solution using Conjugate Gradient */

  init_laplace_matrix(M, sigma, s, L);

  fill_source(b, 2.2, 0.5, L);

  conj_gradient_algorithm(M, f, b, r_hat_targ, L, &num_iter);

  inverse_laplace_operator(check_sol, b, sigma, L, L);

  /* Print Solution and Check Solution */
  std::cout <<"\n\t Conjugate Gradient " << "\t" << std::endl;
  std::cout <<"\n\t Solution: \t\t Check: " << "\t" << std::endl;

  for(i = 0; i < L; i++){
    std::cout << "\t" << f[i] << "\t\t"
    << check_sol[i] << "\t"
    << std::endl;
  }

  /* Solution Conjugate Gradient: Efficient Matrix Implementation  */

  efficient_conj_grad_alg(f, b, sigma, s, r_hat_targ, L, &num_iter);

  /* Print Solution and Check Solution */
  std::cout <<"\n\t Conjugate Gradient: Efficient Matrix Implementation " << "\t\n" << std::endl;
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


  std::cout <<"\n\t Performances Conj Gradient " << "\t\n" << std::endl;

  standard_implemetation.open("../data/standard_timing.dat");


  for(L = L_in; L < L_fin; L += 50){

    M = new double[L * L];
    f = new double[L];
    b = new double[L];

    init_laplace_matrix(M, sigma, s, L);

    fill_source(b, 2.2, 0.5, L);

    t_start = seconds();

    for(i = 0; i < num_rep; ++i)
      conj_gradient_algorithm(M, f, b, r_hat_targ, L, &num_iter);

    t_end = seconds();

    standard_implemetation << L << "\t"
       << t_end - t_start << "\t"
       << std::endl << std::endl;

    std::cout << "\t" << L << "\t"
        << t_end - t_start << "\t"
        << std::endl << std::endl;

    delete [] M;
    delete [] f;
    delete [] b;
}

  standard_implemetation.close();

  std::cout <<"\n\t Performances Efficient Conj Gradient " << "\t\n" << std::endl;

  efficient_implemetation.open("../data/efficient_timing.dat");

  for(L = L_in; L < L_fin; L += 50){

    f = new double[L];
    b = new double[L];

    fill_source(b, 2.2, 0.5, L);

    t_start = seconds();

    for(i = 0; i < num_rep; ++i)
      efficient_conj_grad_alg(f, b, sigma, s, r_hat_targ, L, &num_iter);

    t_end = seconds();

    efficient_implemetation << L << "\t"
       << t_end - t_start << "\t"
       << std::endl << std::endl;

    std::cout << "\t" << L << "\t"
          << t_end - t_start << "\t"
          << std::endl << std::endl;

    delete [] f;
    delete [] b;
}

  efficient_implemetation.close();

}
