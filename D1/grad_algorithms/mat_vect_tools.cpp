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
double * mat_vec_prod(double * A, double * x, int N){

  double * mvprod = new double [N];
  int i, j;

  for (i = 0; i < N; i++)
    mvprod[i] = 0.0;

  for(i = 0; i < N; i++)
    for(j = 0; j < N; j++)
      mvprod[i] += A[i * N + j] * x[j];

  return mvprod;
}
