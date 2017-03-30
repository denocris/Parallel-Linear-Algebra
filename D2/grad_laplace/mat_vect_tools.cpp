#include "mat_vect_tools.h"



double scalar_prod(double* x, double* y, int N){

  double scprod = 0.0;
  int i;

  for(i = 0; i < N; i++)
    scprod += x[i] * y[i];

  return scprod;
}

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



void efficient_mat_vect_prod(double* x, double* row, double sigma, double s, int N){

  // First and last rows
  row[0] = (sigma + 1.) * x[0] + s * x[1] + s * x[N - 1];
  row[N - 1] = s * x[0] + s * x[N - 2] + (sigma + 1.) * x[N - 1];

  // Other rows
  for(int i = 1; i < N - 1; i++){
    row[i] = s * x[i - 1] + (sigma + 1.) * x[i] + s * x[i + 1];
  }
}
