#!/bin/bash
#PBS -l nodes=1:ppn=20
#PBS -l walltime=01:00:00

cd P2.6_seed/D3/
module load openmpi

for ((size=1280000;size<=2560000;size*=2))
do
	for i in  2 4 8 16
	do
		mpirun -np $i ./main_benchmark.x $size  #>> $outfile
	done
done
