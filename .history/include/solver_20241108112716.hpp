#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include "mesh.hpp"
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>
#include "csimple_timer.hpp"
#include<jacobi.hpp>
template <typename T>
/**
 * This calss solves the Laplace equation in 2 D using Jacobi iteration technique.
 */
class Solver
{
public:
    /**
     * Constructor of the solver class
     */
    Solver() {};
    /**
     * Serial solver method.
     * @param Mesh<T>& m: grid instance.
     * @param const size_t &max_steps: Number of iterations.
     * @param const size_t & print_interval: interval for a snapshot of the grid. used to generating the GIF.
     */
    void jacobi(Mesh<T> &m, const size_t &max_steps, const size_t &print_interval);
    /**
     * Parallel solver method.
     * @param Mesh<T>& m: grid instance.
     * @param const size_t &max_steps: Number of iterations.
     * @param const size_t & print_interval: interval for a snapshot of the grid. used to generating the GIF.
     */
    void jacobi_parallel(Mesh<T> &m, const size_t &max_steps, const size_t &print_interval);

private:
    /**
     * Print results to a files
     */
    void print_results(Mesh<T> &m, size_t &counter);
};


template <typename T>
void Solver<T>::print_results(Mesh<T> &m, size_t &counter)
{
    const int width = 3;
    std::ostringstream file_name;
    std::string base_dir = "../data/file_";
    file_name << base_dir << std::setw(width) << std::setfill('0') << counter << ".dat";
    counter++;
    std::ofstream file(file_name.str());
    if (file.is_open())
    {
        file << m << std::endl;
        file.close(); // Close the file
    }
    else
    {
        std::cerr << "Error creating file: " << file_name.str() << std::endl;
    }
}