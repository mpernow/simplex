/* 
Header file for the simplex library.
*/

#ifndef SIMPLEX_H_
#define SIMPLEX_H_

// Definitions
// -----------------------------------

// Parameter for reflection about centroid
#define ALPHA 1.0
// Parameter for further expansion (greater than 1)
#define GAMMA 1.5
// Parameter for contraction
#define RHO 0.5
// Parameter for shrinking
#define SIGMA 0.5

// Cost function
// ----------------------------------

typedef double (*func)(int, double*);


// Nelder--Mead algorithm
// ----------------------------------
int run_simplex(int, double*, double*, func, int);



// Utility functions
// --------------------------------

double **init_simplex(int, double*);

void xisq(func, int, double**, double*);

void update_simplex(int, double*, int, int, int, double**, func);

int max_index(double*, int);

int min_index(double*, int);

void compare(int*, int*, int*, double*, int);


#endif // SIMPLEX_H_
