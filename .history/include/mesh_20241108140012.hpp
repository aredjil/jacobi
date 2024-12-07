#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include "mesh_print.tpp"
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
    Mesh()
    {
        // Default constructor implementation
    }
    /**
     * Class constructor to intialize boundary conditions.
     */
    template <typename F>
    Mesh(int n, F &boundary_conditions, const MPI_Comm &comm);

    /* Print to a file in parallel.
     */
    void print_to_file(int count, std::string filename = "output");
    /* Print local data in a matrix format*/
    std::ostream &print(std::ostream &os, const std::vector<T> &matrix) const;
};
/**
 * Constructor definition to intilize the mesh using the boundary conditions using the function.
 */
template <typename T>
template <typename F>
Mesh<T>::Mesh(int n, F &boundary_conditions, const MPI_Comm &comm) : size(n), comm(comm)
{
    /**
     * Setting the local variables to each processes.
     */
    MPI_Comm_size(comm, &this->npes);
    MPI_Comm_rank(comm, &this->rank);

    this->n_loc = size / this->npes;
    this->rest = size % this->npes;

    if (this->rank < this->rest)
    {
        this->n_loc += 1;
        this->offset = 0;
    }
    else
    {
        this->offset = this->rest;
    }
    /**
     * Resize current field.
     * Here the (+2) is to account for added layer to the bulk
     * of the field defined.
     */
    current_field.resize((this->n_loc + 2) * (this->size + 2), 0.5);
    /**
     * Resize new field.
     */
    new_field.resize((this->n_loc + 2) * (this->size + 2), 0.5);
    boundary_conditions(this->current_field, comm, this->n_loc + 2, this->size + 2, this->rank, this->offset, this->npes);
    new_field = current_field;
};


/*Helper function to print the results either to the screen or to a file*/
template <typename T>
std::ostream &Mesh<T>::print(std::ostream &os, const std::vector<T> &matrix) const
{
    os << std::fixed << std::setprecision(3); // Set precision for all output
    // Iterate over the local grid with proper indexing
    for (int i = 1; i < this->n_loc + 1; i++)
    {
        for (int j = 1; j < this->size + 1; j++)
        {
            os << "\t" << matrix[i * (this->size + 2) + j] << "\t"; // Accessing the matrix element
        }
        os << std::endl; // New line after each row
    }
    return os;
}