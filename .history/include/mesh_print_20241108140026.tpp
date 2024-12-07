#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include"mesh.hpp"
/*Parallel printing to a file.
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
