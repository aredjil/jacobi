#ifndef BOUNDARY_CONDITIONS_H
#define BOUNDARY_CONDITIONS_H
#include <omp.h>
#ifdef _OPENACC
#include<openacc.h>
#endif 
template <typename T>
void boundary_conditions(std::vector<T> &field, const MPI_Comm &comm, const int &n_loc, const int &n, const int &rank, const int &offset, const int &npes)
{
    /**
     * This function applied the boundary conditions on a vector.
     */
    T max_val{100.0};                // maximum value of the oundary conditions
    T step_size = max_val / (n - 1); // Step size
    /**
     * Setting the upper and the left boundaries.
     */
    for (int i = 0; i < n_loc; ++i)
    {
        int i_g = i + ((n_loc - 2) * rank) + offset;
        field[i * n + (n - 1)] = 0.0;
        field[i * n] = max_val - i_g * step_size;
    }
    if (npes == 1)
    {
        for (int i = 0; i < n; ++i) // iterate columns
        {
            field[i + n * (n_loc - 1)] = 0.0;
        }
    }
    if (!rank)
    {
        for (int i = 0; i < n; ++i) // iterate columns
        {
            // Upper boundary.
            field[i] = max_val - (i)*step_size;
        }
    }
    else if (rank == npes - 1)
    {
        for (int i = 0; i < n; ++i) // iterate columns
        {
            field[i + n * (n_loc - 1)] = 0.0;
        }
    }
    
};
#endif // BOUNDARY_CONDITIONS_H