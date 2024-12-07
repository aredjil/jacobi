#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <mpi.h>
#include <iomanip>
#include"mesh.hpp"
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