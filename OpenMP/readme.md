# COMANDOS PARA EJECUTAR EN LOCAL Y EN GUANE



## COMPILAR 
Para compilar el c�digo se hace mediante gcc, se hace de la misma forma tanto para local como para Guane. A continuacion, el comando: 

- gcc -fopenmp omp_fd1d_advection_diffusion_steady.c -o omp_fd1d_advection_diffusion_steady -lm



## EJECUTAR
Para ejecutar en guane y de forma local, se usa el siguiente comando:

- ./omp_fd1d_advection_diffusion_steady


O se puede usar el siguiente comando, que es preferible, el cual nos dice el tiempo real que tomo en ejecutarse. 

- time ./omp_fd1d_advection_diffusion_steady 



## SBATCH
Se usa el siguiente comando para ejecutar mediante SLURM y su archivo .out se llama output_sbatch_fd1d_advection_diffusion_steady.out

- sbatch omp_sbatch_fd1d_advection_diffusion_steady.sbatch



## SALIDA
Los resultados, que obtiene el codigo en omp se imprimiran en el archivo: output_omp_fd1d_advection_diffusion_steady.txt
