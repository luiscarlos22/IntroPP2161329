# COMANDOS PARA EJECUTAR EN GUANE


## MODULO MPI 
Antes de compilar el codigo se debe cargar el modulo mpi

- module load devtools/mpi/openmpi/4.0.1



## COMPILAR 
Para compilar el codigo se hace mediante gcc, se hace de la misma forma tanto para local como para Guane. A continuacion, el comando: 

- mpicc mpi_fd1d_advection_diffusion_steady.c -o mpi_fd1d_advection_diffusion_steady -lm



## EJECUTAR
Para ejecutar en guane y de forma local, se usa el siguiente comando:

- mpirun -np 2 mpi_fd1d_advection_diffusion_steady


O se puede usar el siguiente comando, que es preferible, el cual nos dice el tiempo real que tomo en ejecutarse. 

- time mpirun -np 2 mpi_fd1d_advection_diffusion_steady



## SBATCH
Se usa el siguiente comando para ejecutar mediante SLURM y su archivo .out se llama sbatch_output_fd1d_advection_diffusion_steady.out

- sbatch sbatch_mpi_fd1d_advection_diffusion_steady.sbatch



## SALIDA
Los resultados, que obtiene el codigo en omp se imprimiran en el archivo: output_mpi_fd1d_advection_diffusion_steady_data.txt