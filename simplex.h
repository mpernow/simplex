/* 
Header file for the simplex library.
*/

#ifndef SIMPLEX_H_
#define SIMPLEX_H_

// Definitions
// -----------------------------------

// Here will go a bunch of parameter definitions of the algorithm

// Cost function
// ----------------------------------




// Nelder--Mead algorithm
// ----------------------------------
int run_simplex(int d);



// Utility functions
// --------------------------------

double **init_simplex(int d);

void xisq(int d, double **points_arr, double *current, double *xi2_arr);

void update_simplex(int d, double *xisq_vec, int best, int second, int worst, double *current, double **points_arr);

int max_index(double *arr, int len);

int min_index(double *arr, int len);

void compare(int *best, int *second, int *worst, double *xi2_arr, int d);


#endif // SIMPLEX_H_
