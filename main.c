#include <stdlib.h>
#include <stdio.h>
#include <math.h>

 #include "simplex.h"

double rosenbrock(int d, double* point)
{
    // Implements the Rosenbrock function with d dimensions.
    double f;
    int i;

    f = 0.0;

    for (i = 0; i < d-1; i++)
    {
        f += 100*pow(point[i+1] - pow(point[i], 2), 2) + pow(1 - point[i], 2);
    }
    return f;
}

int main(int argc, char* argv[])
{
    int num_iter = 1000;
    // Takes a starting point as command-line argument
    int d = argc - 1;
    double* start_point = malloc(d * sizeof(double));
    double* final_point = malloc(d * sizeof(double));
    if (d > 0) // if cmd line arguments given
    {
        for (int i = 0; i < d; i++)
        {
            start_point[i] = atof(argv[i + 1]);
        }
    }
    else // no cmd line arguments given
    {
        d = 3;
        start_point[0] = 4;
        start_point[1] = 5;
        start_point[2] = 6;
    }
    printf("%f\n", rosenbrock(d, start_point));
    
    func rosen = &rosenbrock;
    run_simplex(d, start_point, final_point, rosen, num_iter);

    for (int i = 0; i < d; i++)
    {
        printf("%f\t", final_point[i]);
    }
    printf("\n");
    printf("%f\n", rosenbrock(d, final_point));

    free(start_point);
    free(final_point);

    return 0;
}
