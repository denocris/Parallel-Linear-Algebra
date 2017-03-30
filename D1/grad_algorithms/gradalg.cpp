#include "gradalg.h"
#include "mat_vect_tools.h"

void gradient_algorithm(double * A, double * x, double * b, double r_hat_targ,
  int N, int * num_iter){

  //double * r;
  double * t; // A*r_(k-1)
  double r_hat_2, alpha;
  int i;

  double * r = new double [N];

  /* First step */
  for(i = 0; i < N; i++){
    x[i] = 0.0;
    r[i] = b[i];
  }

  * num_iter = 0;

  r_hat_2 = scalar_prod(r, r, N) * (1./scalar_prod(b, b, N));

  while( r_hat_2 > r_hat_targ * r_hat_targ){

    t = mat_vec_prod(A, r, N);
    alpha = scalar_prod(r, r, N) * (1./scalar_prod(r, t, N));

    for(i = 0; i < N; i++){
      x[i] = x[i] + alpha * r[i];
      r[i] = r[i] - alpha * t[i];
    }

    r_hat_2 = scalar_prod(r, r, N) * (1./scalar_prod(b, b, N));

    * num_iter = * num_iter + 1;
    delete [] t;
  }
delete [] r;
}

void conj_gradient_algorithm(double * A, double * x, double * b, double r_hat_targ,
  int N, int * num_iter){

  double * t; // A*p_(k-1)
  double r_hat_2, alpha, beta, scalar_prod_rk;
  int i;

  double * r = new double [N];
  double * p = new double [N];

  /* First step */
  for(i = 0; i < N; i++){
    x[i] = 0.0;
    r[i] = b[i];
    p[i] = b[i];
  }

  * num_iter = 0;

  r_hat_2 = scalar_prod(r, r, N) * (1./ scalar_prod(b, b, N));

  while( r_hat_2 > r_hat_targ * r_hat_targ){

    t = mat_vec_prod(A, p, N);
    alpha = scalar_prod(r, r, N) * (1./ scalar_prod(p, t, N));

    scalar_prod_rk = 1./scalar_prod(r, r, N);

    for(i = 0; i < N; i++){
      x[i] = x[i] + alpha * p[i];
      r[i] = r[i] - alpha * t[i];
    }

    beta = scalar_prod(r, r, N) * scalar_prod_rk;

    for(i = 0; i < N; i++){
      p[i] = r[i] + beta * p[i];
    }

    r_hat_2 = scalar_prod(r, r, N) * (1./scalar_prod(b, b, N));

    * num_iter = * num_iter + 1;
    delete [] t;
  }
delete [] r;
delete [] p;
}
