#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

#include "../header/gradalg.h"
#include "../header/tool_laplace.h"
#include "../header/random_gen.hpp"
#include "../header/inverse_laplace_operator.hpp"

/* MPI for Distributed Memory Parallelization */
#include <mpi.h>

int main(int argc, char ** argv){

  double sigma = 0.6,  s = -0.5;
  double r_hat_targ = 1.e-16;
  double t_start, t_end;

  int i, rk, num_iter, L;

  double * f;
  double * b;

  /* File declaration */

  std::ofstream time_data;

  /* MPI variables declaration */

  int size; // number of processes
  int rank; // name of the process
  int mytag = 16;
  int vect_size, rest, L_local, l_loc_tmp;
  double * results_recv, * b_to_send;
  int * displ, * recv;

  L = atoi(argv[1]);

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);


  /* For Weak Scaling */
  //L = L * size;

  vect_size = L;
  L_local = L / size;
  rest = vect_size % size;

  /* Distribute the rest following the order: rank 0,1,2,... */
  if(rest != 0 && rank < rest)
    L_local++;

  if(rank == 0){
      f = new double[L_local];
      b = new double[L_local];

      displ  = new int[size];
      recv   = new int[size];
      b_to_send = new double[vect_size];

      recv[0] = L_local;
      displ[0] = 0;

      /* Rank 0 fills its local b */
      fill_source(b, 2.2, 0.5, L_local);

      l_loc_tmp = L_local;

      /* Rank 0 sends to the others */
      for(rk = 1; rk < size; rk++){

        if(rest != 0 && rk == rest)
          l_loc_tmp--;

        fill_source(b_to_send, 2.2, 0.5, l_loc_tmp);
        MPI_Send(b_to_send, l_loc_tmp, MPI_DOUBLE, rk, mytag, MPI_COMM_WORLD);

        recv[rk] = l_loc_tmp;
        displ[rk] = displ[rk - 1] + recv[rk - 1];
      }
/* initialization of results_recv pointer need from process 0 to gather the results. */
    results_recv = new double[vect_size];
  }
  else{

    f = new double[L_local];
    b = new double[L_local];

    /* All other processes different from 0 receives in b */
    MPI_Recv(b, L_local, MPI_DOUBLE, 0, mytag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  /* Efficient MPI Implementation of Conj Gradient */

  t_start = MPI_Wtime();

  efficient_conj_grad_alg_MPI(f, b, sigma, s, r_hat_targ, L_local, &num_iter, rank, size);

  t_end = MPI_Wtime();


  /* Takes local f and b and sends them to results_recv and b_to_send in precise positions */
  /* recv: Integer array (of length group size) containing
      the number of elements that are received from each process (significant only at root).*/

  /* displ: Integer array (of length group size). Entry i specifies the displacement relative
      to recvbuf at which to place the incoming data from process i (significant only at root).*/

  MPI_Gatherv(f, L_local, MPI_DOUBLE, results_recv, recv, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gatherv(b, L_local, MPI_DOUBLE, b_to_send, recv, displ, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if(rank == 0){
      /* WEAK SCALING */
/*
      time_data.open("../data/weak_scaling.dat", std::ios_base::app);

      time_data << size << "\t"
         << atoi(argv[1]) << "\t"
         << vect_size << "\t"
         << t_end - t_start << "\t"
         << std::endl;

      time_data.close();
*/
      /* STRONG SCALING */

      time_data.open("../data/strong_scaling.dat", std::ios_base::app);

      time_data << size << "\t"
         << t_end - t_start << "\t"
         << std::endl;

      time_data.close();


    //std::cout <<  "\n\t Benchmark Finished \t\n" << std::endl;

    delete [] results_recv;
    delete [] displ;
    delete [] recv;
    delete [] b_to_send;
  }

  delete [] f;
  delete [] b;
  MPI_Finalize();

  return 0;
}
