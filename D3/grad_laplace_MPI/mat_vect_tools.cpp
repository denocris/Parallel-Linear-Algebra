#include "mat_vect_tools.h"

#include <mpi.h>

double scalar_prod(double * x, double * y, int N){

  double scprod = 0.0;
  double total_scprod = 0.0;
  int i;

  for(i = 0; i < N; i++)
    scprod += x[i] * y[i];

  MPI_Allreduce(&scprod, &total_scprod, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  scprod = total_scprod;

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
  double total_norm;

  for(int i = 0; i < N; i++)
    norm += x[i] * x[i];

  MPI_Allreduce(&norm, &total_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  norm = total_norm;

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

void efficient_mat_vect_prod_MPI(double* x, double* row, double sigma, double s, int N,
  int rank, int size){

  int i;

  int mytag = 16;
  double x_tmp;
  MPI_Request MyReq;

  /* Product of the diagonal elements*/
  for(i = 0; i < N; i++)
    row[i] = (sigma + 1.) * x[i];

  // Sends and receives a message

  /*MPI_Sendrecv(*sendbuf, int sendcount, sendtype, int dest, int sendtag,
              *recvbuf, int recvcount, recvtype,int source, int recvtag,
                    MPI_Comm comm, MPI_Status *status)*/

  MPI_Sendrecv(x, 1, MPI_DOUBLE, (rank+size-1)%size, mytag,
               &x_tmp, 1, MPI_DOUBLE, (rank+1)%size, mytag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  x[0] = x_tmp;

  /* Updating the product */
  for(i = 0; i < N; i++)
    row[i] += s * x[(i+1)%N];

  /* "arrow inversion" */
  MPI_Sendrecv(x, 1, MPI_DOUBLE, (rank+1)%size, mytag,
               &x_tmp, 1, MPI_DOUBLE, (rank+size-1)%size, mytag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  x[0] = x_tmp;

  /* "down shifting" */
  MPI_Sendrecv(&(x[N-1]), 1, MPI_DOUBLE, (rank+1)%size, mytag,
               &x_tmp, 1, MPI_DOUBLE, (rank+size-1)%size, mytag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  x[N-1] = x_tmp;

  for(i = 0; i < N; i++)
    row[i] += s * x[(N+i-1)%N];

  /* "arrow inversion" */
  MPI_Sendrecv(&(x[N-1]), 1, MPI_DOUBLE, (rank+size-1)%size, mytag,
               &x_tmp, 1, MPI_DOUBLE, (rank+1)%size, mytag,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  x[N-1] = x_tmp;
}
