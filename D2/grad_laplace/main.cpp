#include "gradalg.h"
#include "mat_vect_tools.h"
#include "random_gen.hpp"
#include <fstream>
#include <iostream>

int main(){

  double * check_sol;
  double r_hat = 1.e-3;
  double cond_numb = 1.e5;

  int N = 2, num_iter, i;
  //int num_rep = 10;

  double * A = new double [N * N];
  double * b = new double [N];
  double * sol = new double [N];
  double err_impl;
  double err_expl;

  //ofstream scaling_iter;
  std::ofstream errors;

/* Filling matrix with the assignment values */

  A[0] = 3.0;
  A[1] = 1.0;
  A[2] = 1.0;
  A[3] = 2.0;

  b[0] = 1.0;
  b[1] = 3.0;

/* ----------------------------------------------------------------- */

  printf("\n\t Gradient Algorithm \n");

  gradient_algorithm(A, sol, b, r_hat, N, &num_iter);

  check_sol = mat_vec_prod(A, sol, N);

  printf("\n\tsol_x = %lg, sol_y = %lg", sol[0], sol[1]);
  printf("\n\n\tCHECH sol:\n\tcheck_x = %lg, check_y = %lg\n", check_sol[0], check_sol[1]);
  printf("\n\tResult obtained in %d iteration\n", num_iter);

  delete [] check_sol;

/* ----------------------------------------------------------------- */

  printf("\n\t Conj - Gradient Algorithm \n");

  conj_gradient_algorithm(A, sol, b, r_hat, N, &num_iter);

  check_sol = mat_vec_prod(A, sol, N);

  printf("\n\tsol_x = %lg, sol_y = %lg", sol[0], sol[1]);
  printf("\n\n\tCHECH sol:\n\tcheck_x = %lg, check_y = %lg\n", check_sol[0], check_sol[1]);
  printf("\n\tResult obtained in %d iteration\n", num_iter);

// Deallocate pointers
  delete [] check_sol;
  delete [] sol;
  delete [] A;
  delete [] b;

/* ----------------------------------------------------------------- */


  //scaling_iter.open("../data/scaling.dat");

  N = 150;
  cond_numb = 1e5;
  r_hat = 1.e-25;

  A = new double [N * N];
  b = new double [N];
  sol = new double [N];

  fill_defpos_symm_matrix(A, cond_numb, N);
  fill_source(b, 2., 0.5, N);

  errors_conj_grad_alg(A, sol, b, r_hat, N);

  delete [] A;
  delete [] b;
  delete [] sol;


    // std::cout << N << "\t"
    //           << num_iter << "\t"
    //           << err_impl << "\t"
    //           << err_expl << "\t"
    //           << std::abs(err_impl - err_expl)
    //           << std::endl << std::endl;

}
