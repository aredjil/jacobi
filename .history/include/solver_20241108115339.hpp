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
    Solver(Mesh<T> &m, const size_t &max_steps, const size_t &print_interval){
        this->m = m;
        this->max_steps=max_steps;
        this->print_interval=print_interval;
    };
    /**
     * Serial solver method.
     * @param Mesh<T>& m: grid instance.
     * @param const size_t &max_steps: Number of iterations.
     * @param const size_t & print_interval: interval for a snapshot of the grid. used to generating the GIF.
     */
    void jacobi();
    /**
     * Parallel solver method.
     * @param Mesh<T>& m: grid instance.
     * @param const size_t &max_steps: Number of iterations.
     * @param const size_t & print_interval: interval for a snapshot of the grid. used to generating the GIF.
     */
    void jacobi_parallel_blocking();
    /**
     * Parallel non
     */
    void jacobi_parallel_nonblocking();

private:
    /**
     * Print results to a files
     */
    void print_results(Mesh<T> &m, size_t &counter);
};
