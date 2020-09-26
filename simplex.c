/*
simplex.c

Nelder-Mead simplex optimisation algorithm implemented

Author: Marcus Pernow

Date: 2018-10-04

*/


#include <stdlib.h>
#include <stdio.h>
#include "RGrunner.h"


double **init_simplex(int d)
{
	// Initialise simplex of d + 1 points in d dimensions:
	// (1 		1 		1 		1 ...
	//  1+length	1		1		1 ...
	//  1		1+length	1		1 ...
	//  1		1		1 + length	1 ...
	//  ...		...		...		....)



	int i, j;
	double length = 1e-2; // start with 1% deviation

	double **points_arr = calloc(d + 1, sizeof(double*));
	for (i = 0; i < d + 1; i++)
	{
		points_arr[i] = calloc(d, sizeof(double));
		for (j = 0; j < d; j++)
		{
			if (i > 0 && i - 1 == j)
			{
				points_arr[i][j] = 1.0 + length;
			}
			else
			{
				points_arr[i][j] = 1.0;
			}
		}
	}

	return points_arr;
}

void xisq(int d, double **points_arr, double *current, double *xi2_arr)
{
	// Calculate the xi2 for each of the d+1 points in the simplex

	int i, j;
	double temp[d];

	for (i = 0; i < d + 1; i++)
	{
		for (j = 0; j < d; j++)
		{
			temp[j] = points_arr[i][j] * current[j];
		}
		xi2_arr[i] = xisquare(temp);
	}
}

void update_simplex(int d, double *xisq_vec, int best, int second, int worst, double *current, double **points_arr)
{
	// Updates the simplex to the next generation
	int i, j;
	double x_av[d], x_w[d], x_r[d], x_new[d], x_e[d], x_c[d], temp[d];

	double xisq_r, xisq_e, xisq_c, xisq_new;
	// Initialise worst and avg point:
	for (i = 0; i < d; i++)
	{
		x_w[i] = points_arr[worst][i];
		x_av[i] = 0.0;
	}
	// Average of all except worst
	for (j = 0; j < d + 1; j++)
	{
		if (j != worst)
		{
			for (i = 0; i < d; i++)
			{
				x_av[i] += points_arr[j][i]/((double)d);
			}
		}
	}

	// Reflect worst point about x_av
	for (i = 0; i < d; i++)
	{
		x_r[i] = 2.0*x_av[i] - x_w[i];
		temp[i] = x_r[i] * current[i];
	}
	xisq_r = xisquare(temp);

	if (xisq_r < xisq_vec[worst])
	{
		if (xisq_r < xisq_vec[second])
		{
			if (xisq_r > xisq_vec[best])
			{
				// Between best and second best
				// Replace worst by it
				for (i = 0; i < d; i++)
				{
					x_new[i] = x_r[i];
					xisq_new = xisq_r;
				}
			}
			else
			{
				// Better than previous best
				// Extend further
				for (i = 0; i < d; i++)
				{
					x_e[i] = 2.0 * x_r[i] - x_av[i];
					temp[i] = x_e[i] * current[i];
				}
				xisq_e = xisquare(temp);
				if (xisq_e < xisq_vec[best])
				{
					// This is now best
					// Make this the new point
					for (i = 0; i < d; i++)
					{
						x_new[i] = x_e[i];
						xisq_new = xisq_e;
					}
				}
				else
				{
					// x_r was better
					for (i = 0; i < d; i++)
					{
						x_new[i] = x_r[i];
						xisq_new = xisq_r;
					}
				}
			}
		}
		else
		{
			// Beetween worst and second best
			// Contract back towards x_av
			for (i = 0; i < d; i++)
			{
				x_new[i] = x_r[i] - 0.5 * (x_r[i] - x_av[i]);
				temp[i] = x_new[i] * current[i];
			}
			xisq_new = xisquare(temp);
		}
		// Do the replacement
		for (i = 0; i < d; i++)
		{
			points_arr[worst][i] = x_new[i];
			xisq_vec[worst] = xisq_new;
		}
	}
	else
	{
		// Reflected point was worse than previous worst
		// Contract inside simplex
		for (i = 0; i < d; i++)
		{
			x_c[i] = x_av[i] - 0.5*(x_av[i] + x_w[i]);
			temp[i] = x_c[i] * current[i];
		}
		xisq_c = xisquare(temp);
		if (xisq_c < xisq_vec[worst])
		{
			// Successful contraction
			for (i = 0; i < d; i ++)
			{
				points_arr[worst][i] = x_c[i];
			}
			xisq_vec[worst] = xisq_c;
		}
		else
		{
			// No improvement. Shrink the simplex towards best
			for (j = 0; j < d + 1; j++)
			{
				if (j != best)
				{
					for (i = 0; i < d; i++)
					{
						points_arr[j][i] = 0.5*(points_arr[j][i] + points_arr[best][i]);
					}
				}
			}
			xisq(d, points_arr, current, xisq_vec);
		}
	}

}

int max_index(double *arr, int len)
{
	// Finds index of miximum in array
	int k = 0;
	int i;
	double max = arr[k];
	
	for (i = 0; i < len; i++)
	{
		if (arr[i] > max)
		{
			max = arr[i];
			k = i;
		}
	}
	return k;
}

int min_index(double *arr, int len)
{
	// Finds index of miximum in array
	int k = 0;
	int i;
	double min = arr[k];
	
	for (i = 0; i < len; i++)
	{
		if (arr[i] < min)
		{
			min = arr[i];
			k = i;
		}
	}
	return k;
}
void compare(int *best, int *second, int *worst, double *xi2_arr, int d)
{
	// Finds best, second best, worst point in simplex
	*best = min_index(xi2_arr, d+1);
	*worst = max_index(xi2_arr, d+1);
	double *tmp = calloc(d + 1, sizeof(double));
	int i;
	for (i = 0; i < d+1; i++)
	{
		if (i == *best)
		{
			tmp[i] = xi2_arr[*worst];
		}
		else
		{
			tmp[i] = xi2_arr[i];
		}
	}
	*second = min_index(tmp, d+1);
	free(tmp);
}



int run_simplex(int d)
{
	double **points_arr;
	int i, j;
	FILE *fp;
	double current[d];
	double tmp;
	double *xi2_arr;
	int best, second, worst;
	xi2_arr = calloc(d + 1, sizeof(double));

	for (j = 0; j < 1000; j++)
	{
		// Read current best point from file "point.dat"
		fp = fopen("point.dat", "r");
		int line = 0;
		while (fscanf(fp, "%lf", &tmp) != EOF)
		{
			current[line] = tmp;
			line += 1;
			printf("%e\n", tmp);
		}
		
		fclose(fp);

		points_arr = init_simplex(d);

		xisq(d, points_arr, current, xi2_arr);
		compare(&best, &second, &worst, xi2_arr, d);
		
		for (i = 0; i < 500; i++)
		{
			update_simplex(d, xi2_arr, best, second, worst, current, points_arr);
			compare(&best, &second, &worst, xi2_arr, d);
			if (i % 10 == 0)
			{
				printf("%.10f\n", xi2_arr[best]);
			}
		}

		fp = fopen("point.dat", "w");
		for (i = 0; i < d; i++)
		{
			fprintf(fp, "%.10e\n", points_arr[best][i] * current[i]);
		}
		fclose(fp);

		for (i = 0; i < d; i++)
		{
			//printf("%e\n", points_arr[best][i] * current[i]);
		}

		for (i = 0; i < d + 1; i++)
		{
			free(points_arr[i]);
		}
		free(points_arr);
	}

	free(xi2_arr);
	return 0;
}
