#include "mat_vect_tools.h"

/* in this way: (x,y) = \sum_i x_i y_i */
double scalar_prod(double* x, double* y, int N){

  double scprod = 0.0;
  int i;

  for(i = 0; i < N; i++)
    scprod += x[i] * y[i];

  return scprod;
}

/* this function return a pointer to the memory area */
/* that contain the result of the operation (array of size N) */
double * mat_vec_prod(const double * A, const double * x, const int N){

  double * mvprod = new double [N];
  int i, j;

  for (i = 0; i < N; i++)
    mvprod[i] = 0.0;

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      mvprod[i] += A[i * N + j] * x[j];

  return mvprod;
}

double vector_norm(const double * x, const int N){

  double norm = 0.0;

  for(int i = 0; i < N; i++)
    norm += x[i] * x[i];

  return norm;
}
