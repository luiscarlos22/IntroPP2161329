#!/bin/bash

#SBATCH --partition=guane_24_cores  
#SBATCH --job-name=MPI      
#SBATCH --nodes=2           
#SBATCH --ntasks=4          
#SBATCH --cpus-per-task=8   
#SBATCH --output=sbatch_output_fd1d_advection_diffusion_steady.out
module load devtools/mpi/openmpi/4.0.1  # Modulo que contiene el software a usar

if [ -n "$SLURM_CPUS_PER_TASK" ]; then
  MPI_threads=$SLURM_CPUS_PER_TASK
else
  MPI_threads=1
fi

export MPI_NUM_THREADS=$omp_threads
 
mpirun -np 2 mpi_fd1d_advection_diffusion_steady
