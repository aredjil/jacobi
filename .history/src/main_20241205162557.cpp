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

int main(int argc, char** argv)
{
    MPI_Init(nullptr, nullptr);
    int rank, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;

    if(argc < 2){
        if(!rank){
            std::cerr<<"Missing argument"<<std::endl;
             std::cerr<<"Basic usage of the code:"<<std::endl;
             std::cerr<<mpirun -np #proc ./excutable_name <size of the mesh>"
             MPI_Abort(comm, 1);
        }
    }
    int n = std::atoi(argv[1]);
    int print_interval = 1;
    int max_steps = 10;
    
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