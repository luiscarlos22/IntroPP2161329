#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[]);
// double *r8vec_linspace_new(int n, double a, double b);
void timestamp();
// double *trisolve(int n, double a[], double b[]);

/******************************************************************************/

int main(int argc, char *argv[])

/******************************************************************************/
/*
  Purpose:

    FD1D_ADVECTION_DIFFUSION_STEADY solves steady advection diffusion equation.

  Discussion:

    The steady advection diffusion equation has the form:

      v ux - k * uxx = 0

    where V (the advection velocity) and K (the diffusivity) are positive
    constants, posed in the region

      a = 0 < x < 1 = b

    with boundary conditions

      u(0) = 0, u(1) = 1.

    The discrete solution is unreliable when dx > 2 * k / v / ( b - a ).

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    03 May 2014

  Author:

    John Burkardt

  MPI Modification:
    
    07 August 2020 by Luis Carlos Jimenez Arciniegas, Universidad Industrial de Santander lucaja999@gmail.com

    To parallelize with MPI the code that solves the constant advection diffusion equation, the first thing that 
    was done was to pass the code that was in the r8vec_linspace_new () and trisolve () functions to main (), 
    adapting the variables of each function with respect to the main () Some cycles were also omitted from the 
    code, for example declaring the vector x and using it to calculate w, everything is done in the same cycle. 
    MPI_Scatter () and MPI_Gather () have been used to work the for cycle, where in this the declaration of vector
    x and the calculation of w is made, also in this vector a3 is calculated, which is used to calculate u. Where u
    is the solution using the finite difference method, and w is the exact solution of the equation. Because vector 
    a3 turns out to be a 3xn matrix, vector a3 has been divided into three vectors, to facilitate the calculation 
    with MPI.

    Para paralelizar con MPI el código que resuelve la ecuación de difusión de advección constante, lo primero que 
    se hizo fue pasar el código que estaba en las funciones r8vec_linspace_new () y trisolve () a main (), adaptando 
    las variables de cada función con respeto to the main () También se omitieron algunos ciclos del código, por 
    ejemplo, declarando el vector x y usándolo para calcular w, todo se hace en el mismo ciclo. MPI_Scatter () y 
    MPI_Gather () se han utilizado para trabajar el ciclo for, donde en este se realiza la declaración del vector 
    x y el cálculo de w, también en este vector se calcula a3, que se utiliza para calcular u. Donde u es la solución 
    usando el método de diferencia finita, y w es la solución exacta de la ecuación. Debido a que el vector a3 resulta
    ser una matriz 3xn, el vector a3 se ha dividido en tres vectores, para facilitar el cálculo con MPI.

*/
{
    double a;
    double b;
    // double *a3;
    // char command_filename[] = "mpi_fd1d_advection_diffusion_steady_commands.txt";
    // FILE *command_unit;
    char data_filename[] = "output_mpi_fd1d_advection_diffusion_steady_data.txt";
    FILE *data_unit;
    double dx;
    int i;
    int j;
    double k;
    int nx;
    double r;
    double *u;
    double v;
    double *w;
    double *x;
    double xmult;
    int n;
    double *x2;
    double *w2;
    double *a0;
    double *a1;
    double *a2;
    double *a0t;
    double *a1t;
    double *a2t;
    double *ut;
    int numeroProcesadores, idProceso;

    /*
      Physical constants.
    */
    v = 1.0;
    k = 0.05;
    /*
      Spatial discretization.
    */
    nx = 101;
    a = 0.0;
    b = 1.0;
    dx = (b - a) / (double)(nx - 1);
    r = v * (b - a) / k;

    /*
      Iniciar MPI
    */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numeroProcesadores);
    MPI_Comm_rank(MPI_COMM_WORLD, &idProceso);
    if (idProceso == 0)
    {
        timestamp();
        printf("\n");
        printf("MPI_FD1D_ADVECTION_DIFFUSION_STEADY:\n");
        printf("  C version\n");
        printf("\n");
        printf("  Solve the 1D steady advection diffusion equation:,\n");
        printf("    v du/dx - k d2u/dx2 = 0\n");
        printf("  with constant, positive velocity V and diffusivity K\n");
        printf("  over the interval:\n");
        printf("    0.0 <= x <= 1.0\n");
        printf("  with boundary conditions:\n");
        printf("    u(0) = 0, u(1) = 1.\n");
        printf("\n");
        printf("  Use finite differences\n");
        printf("   d u/dx  = (u(t,x+dx)-u(t,x-dx))/2/dx\n");
        printf("   d2u/dx2 = (u(x+dx)-2u(x)+u(x-dx))/dx^2\n");
        /*
        Physical constants.
      */

        printf("\n");
        printf("  Diffusivity K = %g\n", k);
        printf("  Velocity V    = %g\n", v);
        /*
        Spatial discretization.
      */
        printf("  Number of nodes NX = %d\n", nx);
        printf("  DX = %g\n", dx);
        printf("  Maximum safe DX is %g\n", 2.0 * k / v / (b - a));
    }
    /*
    The exact solution to the differential equation is known.
  */

    if ((nx % numeroProcesadores) == 0)
    {
        n = (nx / numeroProcesadores);
    }
    else
    {
        n = (nx / numeroProcesadores) + 1;
    }

    w = (double *)malloc(nx * sizeof(double));
    w2 = (double *)malloc(nx * sizeof(double));
    x = (double *)malloc(nx * sizeof(double));
    x2 = (double *)malloc(nx * sizeof(double));
    a0 = (double *)malloc(nx * sizeof(double));
    a1 = (double *)malloc(nx * sizeof(double));
    a2 = (double *)malloc(nx * sizeof(double));
    u = (double *)malloc(nx * sizeof(double));
    a0t = (double *)malloc(nx * sizeof(double));
    a1t = (double *)malloc(nx * sizeof(double));
    a2t = (double *)malloc(nx * sizeof(double));
    ut = (double *)malloc(nx * sizeof(double));

  /*
    MPI_Scatter is used to divide the vectors into subvectors for each process to use.

    MPI_Scatter se usa para dividir los vectores en subvectores para que cada proceso los use.
  */

    MPI_Scatter(x,       // Matriz que se comparte
        n,               // Numero de columnas que se comparte
        MPI_DOUBLE,      // Tipo de dato a enviar
        x2,              // Vector donde se reciben los datos
        n,               // Numero de columnas que se reciben
        MPI_DOUBLE,      // Tipo de dato a recibir
        0,               // Proceso raiz que envia los datos
        MPI_COMM_WORLD); // Comunicador utilizado
    MPI_Scatter(w, n, MPI_DOUBLE, w2, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(u, n, MPI_DOUBLE, ut, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(a0, n, MPI_DOUBLE, a0t, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(a1, n, MPI_DOUBLE, a1t, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(a2, n, MPI_DOUBLE, a2t, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < n; i++)
    {
        if ((i + (idProceso * n)) < nx)
        {
            x2[i] = ((double)(nx - 1 - i) * a + (double)(i + (idProceso * n)) * b) / (double)(nx - 1);
            w2[i] = (1.0 - exp(r * x2[i])) / (1.0 - exp(r));
            a0t[i] = -v / dx / 2.0 - k / dx / dx;
            a1t[i] = +2.0 * k / dx / dx;
            a2t[i] = +v / dx / 2.0 - k / dx / dx;
            ut[i] = 0.0;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  
  /*
    MPI_Gather is used to get the data from the subvectors and put them together into a single vector.
    
    MPI_Gather se usa para obtener los datos de los subvectores y juntarlos en un solo vector.
  */

    MPI_Gather(x2,              // Dato que envia cada proceso
        n,               // Numero de elementos que se envian
        MPI_DOUBLE,      // Tipo del dato que se envia
        x,               // Vector en el que se recolectan los datos
        n,               // Numero de datos que se esperan recibir por cada proceso
        MPI_DOUBLE,      // Tipo del dato que se recibira
        0,               // proceso que va a recibir los datos
        MPI_COMM_WORLD); // Canal de comunicacion (Comunicador Global)
    MPI_Gather(w2, n, MPI_DOUBLE, w, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(ut, n, MPI_DOUBLE, u, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(a0t, n, MPI_DOUBLE, a0, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(a1t, n, MPI_DOUBLE, a1, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(a2t, n, MPI_DOUBLE, a2, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);

    if (idProceso == 0)
    {
        /*
        Set up the tridiagonal linear system corresponding to the boundary
        conditions and advection-diffusion equation.
      */

        a0[0] = 0.0;
        a1[0] = 1.0;
        a2[0] = 0.0;
        a0[nx - 1] = 0.0;
        a1[nx - 1] = 1.0;
        a2[nx - 1] = 0.0;
        u[0] = 0.0;
       
        u[nx - 1] = 1.0;
        for (i = 0; i < nx; i++)
        {
            if (a1[i] == 0.0)
            {
                fprintf(stderr, "\n");
                fprintf(stderr, "TRISOLVE - Fatal error!\n");
                fprintf(stderr, "  A(%d,2) = 0.\n", i);
                exit(1);
            }
            else if (i > 0)
            {
                xmult = a0[i] / a1[i - 1];
                a1[i] = a1[i] - xmult * a2[i - 1];
                u[i] = u[i] - xmult * u[i - 1];
            }
        }
        u[nx - 1] = u[nx - 1] / a1[nx - 1];
        for (i = nx - 2; 0 <= i; i--)
        {
            u[i] = (u[i] - a2[i] * u[i + 1]) / a1[i];
        }

        /*
        Write data file.
      */
        data_unit = fopen(data_filename, "wt");
        for (j = 0; j < nx; j++)
        {
            fprintf(data_unit, "%g  %g  %g\n", x[j], u[j], w[j]);
        }
        fclose(data_unit);

        printf("\n");
        printf("  Gnuplot data written to file '%s'.\n", data_filename);
        // /*
        //   Write command file.
        // */
        //   command_unit = fopen ( command_filename, "wt" );

        //   fprintf ( command_unit, "set term png\n" );
        //   fprintf ( command_unit, "set output 'fd1d_advection_diffusion_steady.png'\n" );
        //   fprintf ( command_unit, "set grid\n" );
        //   fprintf ( command_unit, "set style data lines\n" );
        //   fprintf ( command_unit, "unset key\n" );
        //   fprintf ( command_unit, "set xlabel '<---X--->'\n" );
        //   fprintf ( command_unit, "set ylabel '<---U(X)--->'\n" );
        //   fprintf ( command_unit, "set title 'Exact: green line, Approx: red dots'\n" );
        //   fprintf ( command_unit, "plot '%s' using 1:2 with points pt 7 ps 2,\\\n", data_filename );
        //   fprintf ( command_unit, "'' using 1:3 with lines lw 3\n" );
        //   fprintf ( command_unit, "quit\n" );

        //   fclose ( command_unit );

        //   printf ( "  Gnuplot commands written to '%s'\n", command_filename );
        /*
        Free memory.
      */
      //free(a3);
        free(a0);
        free(a1);
        free(a2);
        free(u);
        free(w);
        free(x);
        /*
        Terminate.
      */
        printf("\n");
        printf("MPI_FD1D_ADVECTION_DIFFUSION_STEADY\n");
        printf("  Normal end of execution.\n");
        printf("\n");
        timestamp();
    }
    MPI_Finalize();
    return 0;
}
/******************************************************************************/

// double *r8vec_linspace_new(int n, double a, double b)

// /******************************************************************************/
// /*
//   Purpose:

//     R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.

//   Discussion:

//     An R8VEC is a vector of R8's.

//     4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.

//     In other words, the interval is divided into N-1 even subintervals,
//     and the endpoints of intervals are used as the points.

//   Licensing:

//     This code is distributed under the GNU LGPL license.

//   Modified:

//     29 March 2011

//   Author:

//     John Burkardt

//   Parameters:

//     Input, int N, the number of entries in the vector.

//     Input, double A, B, the first and last entries.

//     Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
// */
// {
//   int i;
//   double *x;

//   x = (double *)malloc(n * sizeof(double));

//   if (n == 1)
//   {
//     x[0] = (a + b) / 2.0;
//   }
//   else
//   {
//     for (i = 0; i < n; i++)
//     {
//       x[i] = ((double)(n - 1 - i) * a + (double)(i)*b) / (double)(n - 1);
//     }
//   }
//   return x;
// }
/******************************************************************************/

void timestamp()

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM

  Licensing:

    This code is distributed under the GNU LGPL license.

  Modified:

    24 September 2003

  Author:

    John Burkardt

  Parameters:

    None
*/
{
    #define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    fprintf(stdout, "%s\n", time_buffer);

    return;
    #undef TIME_SIZE
}
/******************************************************************************/

// double *trisolve(int n, double a[], double b[])

// /******************************************************************************/
// /*
//   Purpose:

//     TRISOLVE factors and solves a tridiagonal system.

//   Discussion:

//     The three nonzero diagonals of the N by N matrix are stored as 3
//     columns of an N by 3 matrix.

//   Example:

//     Here is how a tridiagonal matrix of order 5 would be stored:

//        *  A11 A12
//       A21 A22 A23
//       A32 A33 A34
//       A43 A44 A45
//       A54 A55  *

//   Licensing:

//     This code is distributed under the GNU LGPL license.

//   Modified:

//     03 May 2014

//   Author:

//     John Burkardt

//   Parameters:

//     Input, int N, the order of the linear system.

//     Input/output, double A[N*3].
//     On input, the tridiagonal matrix.
//     On output, the data in these vectors has been overwritten
//     by factorization information.

//     Input, double B[N], the right hand side of the linear system.

//     Output, double TRISOLVE[N], the solution of the linear system.
// */
// {
//   int i;
//   double *x;
//   double xmult;
//   /*
//   The diagonal entries can't be zero.
// */
//   for (i = 0; i < n; i++)
//   {
//     if (a[i + 1 * n] == 0.0)
//     {
//       fprintf(stderr, "\n");
//       fprintf(stderr, "TRISOLVE - Fatal error!\n");
//       fprintf(stderr, "  A(%d,2) = 0.\n", i);
//       exit(1);
//     }
//   }

//   x = (double *)malloc(n * sizeof(double));

//   for (i = 0; i < n; i++)
//   {
//     x[i] = b[i];
//   }

//   for (i = 1; i < n; i++)
//   {
//     xmult = a[i + 0 * n] / a[i - 1 + 1 * n];
//     a[i + 1 * n] = a[i + 1 * n] - xmult * a[i - 1 + 2 * n];
//     x[i] = x[i] - xmult * x[i - 1];
//   }

//   x[n - 1] = x[n - 1] / a[n - 1 + 1 * n];
//   for (i = n - 2; 0 <= i; i--)
//   {
//     x[i] = (x[i] - a[i + 2 * n] * x[i + 1]) / a[i + 1 * n];
//   }

//   return x;
// }
