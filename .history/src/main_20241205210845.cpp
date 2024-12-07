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
#include <omp.h>
#include "boundary_condition.hpp"
#include "csimple_timer.hpp"
#include <cstdlib>
#ifdef _OPENACC
#include <openacc.h>
#endif
int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cerr << "Missing argument" << std::endl;
        std::cerr << "Basic usage of the code:" << std::endl;
        std::cerr << "mpirun -np #proc ./excutable_name <size of the mesh>" << std::endl;
        exit(1);
    }
    MPI_Init(nullptr, nullptr);
    int rank, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    // const char* omp_env = std::getenv("OMP_NUM_THREADS");
    // if (omp_env)
    //     std::cout << "OMP_NUM_THREADS: " << omp_env << std::endl;
    // else
    //     std::cout << "OMP_NUM_THREADS is not set." << std::endl;

    const int n = std::atoi(argv[1]);
    const int print_interval = 1;
    const int max_steps = 10;

    Mesh<double> m(n, boundary_conditions<double>, comm);

    Solver<double> solver(m, max_steps, print_interval);
#ifdef BLOCKING
    std::cout<<"Running serial solver"<<std::endl;
    solver.jacobi_parallel_blocking();
#elif defined(NON_BLOCKING)
    std::cout<<"Running parallel solver with non-blocking"<<std::endl;
    solver.jacobi_parallel_nonblocking();
#elif defined(HYBRID)
    std::cout << "Running hybrid solver" << std::endl;
    solver.jacobi_hybrid();
#elif defined(SERIAL)
    std::cout << "Running serial solver" << std::endl;
    solver.jacobi();
#else
    std::cout << "There have been an issue with the directives" << std::endl;
#endif
    if (!rank)
    {
        SimpleTimer::print_timing_results();
    }
    MPI_Finalize();
    return 0;
}