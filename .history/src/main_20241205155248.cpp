#include "mesh.hpp"
#include "jacobi.tpp"
#include "solver.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <mpi.h>
#include "boundary_condition.hpp"
#include "csimple_timer.hpp"

int main()
{
    MPI_Init(nullptr, nullptr);

    if(argc < 2){
        
    }
    int n = std::atoi(argv[1]);
    int print_interval = 1;
    int max_steps = 10;
    int rank, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    MPI_Comm comm;
    comm = MPI_COMM_WORLD;

    
    Mesh<double> m(n, boundary_conditions<double>, comm);

    // std::cout << m << std::endl;
    Solver<double> solver(m, max_steps, print_interval);
    #ifdef BLOCKING 
    solver.jacobi_parallel_blocking();
    #elif NON_BLOCKING
    solver.jacobi_parallel_nonblocking();
    #elif HYBRID 
    solver.jacobi_hybrid();
    #else
    solver.jacobi();
    #endif 
    solver.jacobi_parallel_blocking();
    if (!rank)
        SimpleTimer::print_timing_results();

    MPI_Finalize();
    return 0;
}