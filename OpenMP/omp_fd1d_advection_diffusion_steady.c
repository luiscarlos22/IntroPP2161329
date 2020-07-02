#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

int main();
//double *r8vec_linspace_new(int n, double a, double b);
void timestamp();
//double *trisolve(int n, double a[], double b[]);

/******************************************************************************/

int main()

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
  
  OpenMP Modification:
    
    01 July 2020 by Luis Carlos Jimenez Arciniegas, Universidad Industrial de Santander lucaja999@gmail.com
    This OpenMP Modification makes a parallelization of the original Code.
    To parallelize the code that solves the constant advection diffusion equation, the first thing that was 
    done was to pass the code that was in the functions r8vec_linspace_new () and trisolve () to main (), 
    adapting the variables of each function with respect to the main () Some cycles were also omitted in the 
    code, for example declaring vector x and using it to compute w, everything is done in the same cycle. 
    Because the calculation of u with respect to w are two totally different processes, and do not depend on 
    their variables, both processes would be carried out in parallel. Where u is the solution using the finite 
    difference method, and w is the exact solution of the equation.

    Para paralelizar el código que resuelve la ecuación de difusión de advección constante, lo primero que se 
    hizo fue pasar el código que estaba en las funciones r8vec_linspace_new () y trisolve () a main (), adaptando 
    las variables de cada función con respecto a la main () Algunos ciclos también se omitieron en el código, 
    por ejemplo, declarando el vector x y usándolo para calcular w, todo se hace en el mismo ciclo. Debido a que 
    el cálculo de u con respecto a w son dos procesos totalmente diferentes, y no dependen de sus variables, 
    ambos procesos se llevarían a cabo en paralelo. Donde u es la solución usando el método de diferencia finita, 
    y w es la solución exacta de la ecuación.
*/
{
  double a;
  double *a3;
  double b;
  // char command_filename[] = "fd1d_advection_diffusion_steady_commands.txt";
  // FILE *command_unit;
  char data_filename[] = "output_omp_fd1d_advection_diffusion_steady.txt";
  FILE *data_unit;
  double dx;
  double xmult;
  int i;
  int j;
  double k;
  int nx;
  double r;
  double *u;
  double v;
  double *w;
  double *x;

  timestamp();
  printf("\n");
  printf("OMP_FD1D_ADVECTION_DIFFUSION_STEADY:\n");
  printf("  C/OpenMP version\n");
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
  v = 1.0;
  k = 0.05;
  printf("\n");
  printf("  Diffusivity K = %g\n", k);
  printf("  Velocity V    = %g\n", v);
  /*
  Spatial discretization.
*/
  nx = 101;
  a = 0.0;
  b = 1.0;
  dx = (b - a) / (double)(nx - 1);
  x = (double *)malloc(nx * sizeof(double));
  printf("  Number of nodes NX = %d\n", nx);
  printf("  DX = %g\n", dx);
  printf("  Maximum safe DX is %g\n", 2.0 * k / v / (b - a));

  /*
  A3 is declared going through the cycle, in parallel.

  A3 se declara atravesando el ciclo, en paralelo.
*/
  a3 = (double *)malloc(nx * 3 * sizeof(double));
  u = (double *)malloc(nx * sizeof(double));
  r = v * (b - a) / k;
  w = (double *)malloc(nx * sizeof(double));

#pragma omp parallel private(i)
  {
    /*
  A3 is declared going through the cycle, in parallel.

  A3 se declara atravesando el ciclo, en paralelo.
*/
#pragma omp for nowait
    for (i = 1; i < nx - 1; i++)
    {

      a3[i + 0 * nx] = -v / dx / 2.0 - k / dx / dx;
      a3[i + 1 * nx] = +2.0 * k / dx / dx;
      a3[i + 2 * nx] = +v / dx / 2.0 - k / dx / dx;
      u[i] = 0.0;
    }

#pragma omp barrier
    /*
    The threads are synchronized, and then a single thread is responsible for obtaining u, 
    while the other parallel threads are responsible for obtaining w.

    Los hilos se sincronizan, y luego un solo hilo es responsable de obtener u, 
    mientras que los otros hilos en paralelo son responsables de obtener w.
*/

#pragma omp single
    {
      a3[0 + 1 * nx] = 1.0;
      u[0] = 0.0;
      a3[nx - 1 + 1 * nx] = 1.0;
      u[nx - 1] = 1.0;
      for (i = 0; i < nx; i++)
      {
        if (a3[i + 1 * nx] == 0.0)
        {
          fprintf(stderr, "\n");
          fprintf(stderr, "TRISOLVE - Fatal error!\n");
          fprintf(stderr, "  A(%d,2) = 0.\n", i);
          exit(1);
        }
        else if (i > 0)
        {
          xmult = a3[i + 0 * nx] / a3[i - 1 + 1 * nx];
          a3[i + 1 * nx] = a3[i + 1 * nx] - xmult * a3[i - 1 + 2 * nx];
          u[i] = u[i] - xmult * u[i - 1];
        }
      }
      u[nx - 1] = u[nx - 1] / a3[nx - 1 + 1 * nx];
      for (i = nx - 2; 0 <= i; i--)
      {
        u[i] = (u[i] - a3[i + 2 * nx] * u[i + 1]) / a3[i + 1 * nx];
      }
    }

    if (nx == 1)
    {
      x[0] = (a + b) / 2.0;
      w[0] = (1.0 - exp(r * x[0])) / (1.0 - exp(r));
    }
    else
    {
#pragma omp for nowait
      for (i = 0; i < nx; i++)
      {
        x[i] = ((double)(nx - 1 - i) * a + (double)(i)*b) / (double)(nx - 1);
        w[i] = (1.0 - exp(r * x[i])) / (1.0 - exp(r));
      }
    }
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
  /*
  Write command file.
*/
  // command_unit = fopen(command_filename, "wt");

  // fprintf(command_unit, "set term png\n");
  // fprintf(command_unit, "set output 'fd1d_advection_diffusion_steady.png'\n");
  // fprintf(command_unit, "set grid\n");
  // fprintf(command_unit, "set style data lines\n");
  // fprintf(command_unit, "unset key\n");
  // fprintf(command_unit, "set xlabel '<---X--->'\n");
  // fprintf(command_unit, "set ylabel '<---U(X)--->'\n");
  // fprintf(command_unit, "set title 'Exact: green line, Approx: red dots'\n");
  // fprintf(command_unit, "plot '%s' using 1:2 with points pt 7 ps 2,\\\n", data_filename);
  // fprintf(command_unit, "'' using 1:3 with lines lw 3\n");
  // fprintf(command_unit, "quit\n");

  // fclose(command_unit);

  // printf("  Gnuplot commands written to '%s'\n", command_filename);
  /*
  Free memory.
*/
  free(a3);
  free(u);
  free(w);
  free(x);
  /*
  Terminate.
*/
  printf("\n");
  printf("OMP_FD1D_ADVECTION_DIFFUSION_STEADY\n");
  printf("  Normal end of execution.\n");
  printf("\n");
  timestamp();
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
