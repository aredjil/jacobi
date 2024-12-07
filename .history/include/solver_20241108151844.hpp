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
    int max_steps; 
    int print_interval;
    /**
     * Constructor of the solver class
     */
    Solver(Mesh<T> &m, const int &max_steps, const int &print_interval){
        this->m = m;
        this->max_steps=max_steps;
        this->print_interval=print_interval;
    };
    /**Serial jacobi solver */
    void jacobi();
    /*parallel blocking jacobi solver*/
    void jacobi_parallel_blocking();
    /*Parallel non-blocking jacobi solver */
    void jacobi_parallel_nonblocking();
    /*Parallel Hybrid blocking jacobi solver*/
    jacobi_parallel_blocking

private:
    /**
     * Print results to a files
     */
    void print_results(Mesh<T> &m, int &counter);
};
