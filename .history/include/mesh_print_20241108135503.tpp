#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include"mesh.hpp"
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

/**
 * Parallel printing to a file.
 */
template <typename T>
void Mesh<T>::print_to_file(int count, std::string filename)
{
    std::string results_dir = "../data/";    // Base directory where to save results
    std::string idx = std::to_string(count); // Convert the count from int to string count starts in Solver class.
    const int width = 3;
    std::ostringstream file_name;
    file_name << results_dir << filename << std::setw(width) << std::setfill('0') << idx << ".dat";
    std::ofstream file(file_name.str());

    if (file.is_open()) // Check if the file was opened correctly.
    {
        if (!this->rank)
        {
            print(file, new_field);
            std::vector<T> new_field((this->n_loc + 2) * (this->size + 2));
            for (int count = 1; count < this->npes; count++)
            {
                MPI_Recv(new_field.data(), (this->n_loc + 2) * (this->size + 2) * sizeof(T), MPI_BYTE, count, count, this->comm, MPI_STATUS_IGNORE);
                print(file, new_field);
            }
        }
        else
        {
            MPI_Send(this->new_field.data(), ((this->n_loc + 2) * (this->size + 2)) * sizeof(T), MPI_BYTE, 0, this->rank, comm);
        }
    }
    else
    {
        std::cerr << "There have been a problem opening the file: " << file_name.str() << std::endl;
    }

    file.close();
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