/**
 * This file contains the class Solver.
 */
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
#include <jacobi.tpp>
template <typename T>
/**
 * This calss solves the Laplace equation in 2 D using Jacobi iteration technique.
 */
class Solver
{
public:
    Mesh<T> m;
    size_t max_steps; 
    size_t print_interval;
    /**
     * Constructor of the solver class
     */
    Solver(Mesh<T> &m){
        this->m = m;
    };
    /**
     * Serial solver method.
     * @param Mesh<T>& m: grid instance.
     * @param const size_t &max_steps: Number of iterations.
     * @param const size_t & print_interval: interval for a snapshot of the grid. used to generating the GIF.
     */
    void jacobi(const size_t &max_steps, const size_t &print_interval);
    /**
     * Parallel solver method.
     * @param Mesh<T>& m: grid instance.
     * @param const size_t &max_steps: Number of iterations.
     * @param const size_t & print_interval: interval for a snapshot of the grid. used to generating the GIF.
     */
    void jacobi_parallel(const size_t &max_steps, const size_t &print_interval);

private:
    /**
     * Print results to a files
     */
    void print_results(Mesh<T> &m, size_t &counter);
};
