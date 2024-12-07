#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include"mesh"
template <typename T>
class Mesh
{
public:
    // MPI communicator to be passed in the constructor
    MPI_Comm comm;
    // Current rank
    int rank;
    // Number of porcesses used
    int npes;
    // Local number of rows in each processes
    int n_loc;
    // The rest of rows if the first dimension of the matrix isn't a multiple of npes.
    int rest;
    // Offset = {0, 1} to be used in printing
    int offset;
    // Size of the field/mesh/matrix the terms are used interchangbly in this code :).
    int size; // size of the mesh.
    // A contigious array to hold the values of the field at t.
    std::vector<T> current_field;
    // A contigious array to hold the values of the field at t+1.
    std::vector<T> new_field;
    /*Default Mesh constructor */
    Mesh() {
        // Default constructor implementation
    }
    /**
     * Class constructor to intialize boundary conditions.
     */
    template <typename F>
    Mesh(int n, F &boundary_conditions, const MPI_Comm &comm);
    
    /**
     * Print to a file in parallel.
     */
    void print_to_file(int count, std::string filename="output");
    /**
     * Print local data in a matrix format
     */
    std::ostream& print(std::ostream &os, const std::vector<T> &matrix) const;
};
