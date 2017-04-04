#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=01:00:00

cd P2.6_seed/D3/grad_laplace_MPI/src_bench
module load openmpi


#for ((size=100;size<=1000;size+=100))
size=500000
do
	for i in  1 2 4 8 16
	do
		mpirun -np $i ./main_benchmark.x $size  #>> $outfile
	done
done
