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


void errors_conj_grad_alg(double * A, double * x, double * b, double r_hat_targ,
  int N){

  double * t;
  double * Ax;
  double r_hat_2, alpha, beta;
  double scalar_prod_rk, square_norm_of_b;
  double e_impl, e_expl;
  int i, num_iter;

  double * err_impl = new double [N];
  double * err_expl = new double [N];
  double * p = new double [N];

  /* First step */
  for(i = 0; i < N; i++){
    x[i] = 0.0;
    err_impl[i] = b[i];
    p[i] = b[i];
  }

  num_iter = 0;

  square_norm_of_b = scalar_prod(b, b, N);

  r_hat_2 = scalar_prod(err_impl, err_impl, N) * (1./ square_norm_of_b);


  std::ofstream errors;
  errors.open("../data/errors.dat");

  while( r_hat_2 > r_hat_targ * r_hat_targ){

    t = mat_vec_prod(A, p, N);
    alpha = scalar_prod(err_impl, err_impl, N) * (1./ scalar_prod(p, t, N));

    scalar_prod_rk = 1./scalar_prod(err_impl, err_impl, N);

    Ax = mat_vec_prod(A, x, N);

    /* Explicit Error */
    for(i = 0; i < N; i++){
      err_expl[i] = b[i] - Ax[i];
    }

    /* Implicit Error */
    for(i = 0; i < N; i++){
      x[i] = x[i] + alpha * p[i];
      err_impl[i] = err_impl[i] - alpha * t[i];
    }

    beta = scalar_prod(err_impl, err_impl, N) * scalar_prod_rk;

    for(i = 0; i < N; i++){
      p[i] = err_impl[i] + beta * p[i];
    }

    r_hat_2 = scalar_prod(err_impl, err_impl, N) * (1./ square_norm_of_b);

    for(i = 0; i < N; i++){
      p[i] = err_impl[i] + beta * p[i];
    }

    /* Error norm computation */


    e_impl = std::pow(vector_norm(err_impl, N) / square_norm_of_b, 0.5);
    e_expl = std::pow(vector_norm(err_expl, N) / square_norm_of_b, 0.5);

    errors << num_iter << "\t"
           << e_impl << "\t"
           << e_expl << "\t"
           << std::abs(e_impl - e_expl)
           << std::endl << std::endl;

    num_iter = num_iter + 1;
    delete [] t;
    delete [] Ax;

  }

errors.close();
delete [] p;
delete [] err_impl;
delete [] err_expl;
}


void efficient_conj_grad_alg(double * x, double * b,
  double sigma, double s, double r_hat_targ, int N, int * num_iter){

  double r_hat_2, alpha, beta;
  double scalar_prod_rk, square_norm_of_b;
  int i;

  double * r = new double [N];
  double * p = new double [N];
  double * t = new double [N];


  /* First step */
  for(i = 0; i < N; i++){
    x[i] = 0.0;
    r[i] = b[i];
    p[i] = b[i];
  }

  * num_iter = 0;

  square_norm_of_b = scalar_prod(b, b, N);

  r_hat_2 = scalar_prod(r, r, N) * (1./ square_norm_of_b);

  while( r_hat_2 > r_hat_targ * r_hat_targ){

    efficient_mat_vect_prod(p, t, sigma, s, N);
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
  }
delete [] t;
delete [] r;
delete [] p;
}
