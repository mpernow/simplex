#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// #include "simplex.h"

double rosenbrock(int d, double* point)
{
    // Implements the Rosenbrock function with d dimensions.
    double f;
    int i;
    for (i = 0; i < d-1; i++)
    {
        f += 100*pow(point[i+1] - pow(point[i], 2), 2) + pow(1 - point[i], 2);
    }
    return f;
}

int main(int argc, char* argv[])
{
    // Takes a starting point as command-line argument
    const int d = argc - 1;
    double* start_point = malloc(d * sizeof(double));
    double* final_point = malloc(d * sizeof(double));
    for (int i = 0; i < d; i++)
    {
        start_point[i] = atof(argv[i + 1]);
    }
    printf("%f\n", rosenbrock(d, start_point));
    
    // Call simplex here!
    final_point = run_simplex(d, start_point, &rosenbrock);

    free(start_point);

    return 0;
}
