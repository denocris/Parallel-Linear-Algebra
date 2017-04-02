#include "gradalg.h"
#include "mat_vect_tools.h"

#include "random_gen.hpp"

int main(){

  double * check_sol;
  double r_hat = 1.e-10;
  double cond_numb = 1.e6;

  int N = 2, num_iter, i;
  int num_rep = 10;

  double * A = new double [N * N];
  double * b = new double [N];
  double * sol = new double [N];

  FILE* scaling_iter;

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

  scaling_iter = fopen("../data/scaling.dat", "w");

  for(i = 0; i < 10; i++){
    for(N = 10; N < 501; N += 20){

      A = new double [N * N];
      b = new double [N];
      sol = new double [N];

      fill_defpos_symm_matrix(A, cond_numb, N);
      fill_source(b, 2., 0.5, N);

      conj_gradient_algorithm(A, sol, b, r_hat, N, &num_iter);

      fprintf(scaling_iter, "%d\t%d\n", N, num_iter);

      delete [] A;
      delete [] b;
      delete [] sol;
    }
  }

  fclose(scaling_iter);
/* ----------------------------------------------------------------- */


int Nfixed = 500;
A = new double [Nfixed * Nfixed];
b = new double [Nfixed];
sol = new double [Nfixed];

scaling_iter = fopen("../data/scaling2.dat", "w");

for(i = 0; i < num_rep; i++){
  for(cond_numb = 1250; cond_numb < 1250 * 6 + 1; cond_numb += 1250){

    fill_defpos_symm_matrix(A, cond_numb, Nfixed);
    fill_source(b, 2., 0.5, Nfixed);

    conj_gradient_algorithm(A, sol, b, r_hat, Nfixed, &num_iter);

    fprintf(  scaling_iter, "%lg\t%d\n", cond_numb, num_iter);

  }
}

delete [] A;
delete [] b;
delete [] sol;

fclose(scaling_iter);


}
