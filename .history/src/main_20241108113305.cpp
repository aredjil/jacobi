#include "mesh.hpp"
#include"jacobi.hpp"
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

int main(int argc, char **argv)
{
    MPI_Init(nullptr, nullptr);
    size_t n = 120;
    int print_interval = 1;
    int max_steps = 10;
    int rank, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    MPI_Comm comm;
    comm = MPI_COMM_WORLD;
    Mesh<double> m(n, boundary_conditions<double>, comm);
    // std::cout << m << std::endl;
    Solver<double> solver;
    solver.jacobi_parallel(m, max_steps, print_interval);
    // Send last local row to the next process and receive into the top halo row from the previous process
    if (!rank)
        SimpleTimer::print_timing_results();

    MPI_Finalize();
    return 0;
}