/**
 * This file contains the methods for the class Solver.
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
#include "solver.hpp"

/**
 * Serial solver.
 */
template <typename T>
void Solver<T>::jacobi()
{
    size_t n = m.size;
    size_t counter = 1;
    for (size_t k = 0; k < max_steps; k++) // Iterate over the total number of steps
    {
        for (size_t i = 1; i < n - 1; i++) // Iterate the rows excluding the boundaries.
        {
            for (size_t j = 1; j < n - 1; j++) // Iterate the columns excluding the boundaries.
            {
                // Update rule.
                m.new_field[i * n + j] = (1.0 / 4) * (m.current_field[(i - 1) * n + j] + m.current_field[(i + 1) * n + j] + m.current_field[i * n + j - 1] + m.current_field[i * n + j + 1]);
            }
        }
        // Swap the current field values with the new field values before the next iteration
        m.new_field.swap(m.current_field);

        if (k % print_interval == 0)
        {
// TODO change this to be inside ifdef
#ifdef PRINT
            print_results(m, counter);
#endif
        }
    }
}
/**
 * Parallel solver.
 */
template <typename T>
void Solver<T>::jacobi_parallel()
{
    size_t n = m.size;
    size_t counter = 1;
    /**
     * Send last row to the next process,
     * and the first row to the previous processes.
     */
    // int rank, npes;
    // MPI_Comm_rank(comm, &rank);
    // MPI_Comm_size(comm, &npes);

    for (size_t k = 0; k < max_steps; k++) // Iterate over the total number of steps
    {

#ifdef B
        // send to the next proc
        MPI_Sendrecv(
            m.current_field.data() + (m.size + 2) * (m.n_loc), (m.size + 2) * sizeof(T), MPI_BYTE,
            m.rank != m.npes - 1 ? m.rank + 1 : MPI_PROC_NULL, 0, m.current_field.data(), (m.size + 2) * sizeof(T),
            MPI_BYTE, !m.rank ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, MPI_STATUS_IGNORE);
        // send to the ghost cell of the previous proc
        MPI_Sendrecv(
            m.current_field.data() + (m.size + 2), (m.size + 2) * sizeof(T), MPI_BYTE,
            !m.rank ? MPI_PROC_NULL : m.rank - 1, 0, m.current_field.data() + (m.size + 2) * (m.n_loc + 1), (m.size + 2) * sizeof(T),
            MPI_BYTE, m.rank != m.npes - 1 ? m.rank + 1 : MPI_PROC_NULL, 0, m.comm, MPI_STATUS_IGNORE);
#elif defined(NONB)
        MPI_Request requests[4]; // Two sends and two receives

        // Send to the next process and receive from the previous process
        MPI_Irecv(m.current_field.data(), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == 0 ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, &requests[0]);
        MPI_Irecv(m.current_field.data() + (m.size + 2) * (m.n_loc + 1), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == m.npes - 1 ? MPI_PROC_NULL : m.rank + 1, 0, m.comm, &requests[2]);
        // Send to the previous process and receive from the next process

        MPI_Isend(m.current_field.data() + (m.size + 2) * m.n_loc, (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == m.npes - 1 ? MPI_PROC_NULL : m.rank + 1, 0, m.comm, &requests[1]);
        MPI_Isend(m.current_field.data() + (m.size + 2), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == 0 ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, &requests[3]);

        // Wait for all non-blocking operations to complete
        // MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
#else
        if (!m.rank)
        {
            std::cout << "There has been an issue" << std::endl;
        }
        break;
#endif
        {
            SimpleTimer t("Commputution Duration");

            for (size_t i = 1; i < m.n_loc + 1; i++) // Iterate the rows excluding the boundaries.
            {
                for (size_t j = 1; j < n + 1; j++) // Iterate the columns excluding the boundaries.
                {
                    // Update rule.
                    m.new_field[i * (n + 2) + j] = (1.0 / 4) * (m.current_field[(i - 1) * (n + 2) + j] + m.current_field[(i + 1) * (n + 2) + j] + m.current_field[i * (n + 2) + j - 1] + m.current_field[i * (n + 2) + j + 1]);
                }
            }
        }
#ifdef NONB
        {
            SimpleTimer t("Communication Duration");
            MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
        }
#endif
        // Swap the current field values with the new field values before the next iteration
        m.new_field.swap(m.current_field);

        if (k % print_interval == 0)
        {
#ifdef PRINT
            // print_results(m, counter);
            m.print_to_file(counter);
            counter++;
#endif
        }
    }
}

// Serial printing function. 
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