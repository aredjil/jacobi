/**
 * This file contains the methods for the class Solver.
 */
#ifndef JACOBI_H
#define JACOBI_H
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
#include "mpi_related_funcs.tpp"
/*Methods definitions for the solver class*/
/* Serial jacobi solver*/
template <typename T>
void Solver<T>::jacobi()
{
    int n = m.size;
    SimpleTimer t("Timing for serial jacobi solver");
    for (int k = 0; k < max_steps; k++) // Iterate over the total number of steps
    {
        for (int i = 1; i < n - 1; ++i) // Iterate the rows excluding the boundaries.
        {
            for (int j = 1; j < n - 1; ++j) // Iterate the columns excluding the boundaries.
            {
                // Update rule.
                m.new_field[i * n + j] = (1.0 / 4) * (m.current_field[(i - 1) * n + j] + m.current_field[(i + 1) * n + j] + m.current_field[i * n + j - 1] + m.current_field[i * n + j + 1]);
            }
        }
        // Swap the current field values with the new field values before the next iteration
        std::swap(m.new_field, m.current_field);

        if (k % print_interval == 0)
        {
#ifdef PRINT
            print_results(m, counter);
#endif
        }
    }
}
/**
 * Parallel blokcing solver.
 */
template <typename T>
void Solver<T>::jacobi_hybrid()
{
    int n = m.size;
    /**
     * Send last row to the next process,
     * and the first row to the previous processes.
     */
    for (int k = 0; k < max_steps; ++k) // Iterate over the total number of steps
    {

        MPI_Request requests[4]; // Two sends and two receives

        // Send to the next process and receive from the previous process
        MPI_Irecv(m.current_field.data(), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == 0 ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, &requests[0]);
        MPI_Irecv(m.current_field.data() + (m.size + 2) * (m.n_loc + 1), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == m.npes - 1 ? MPI_PROC_NULL : m.rank + 1, 0, m.comm, &requests[1]);
        // Send to the previous process and receive from the next process

        MPI_Isend(m.current_field.data() + (m.size + 2) * m.n_loc, (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == m.npes - 1 ? MPI_PROC_NULL : m.rank + 1, 0, m.comm, &requests[2]);
        MPI_Isend(m.current_field.data() + (m.size + 2), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == 0 ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, &requests[3]);
        {
            SimpleTimer t("Computation timing for serial jacobi solver");
            #pragma omp parallel for collapse(2)
            for (int i = 2; i < m.n_loc; i++) // Iterate the rows excluding the boundaries.
            {
                for (int j = 1; j < n + 1; ++j) // Iterate the columns excluding the boundaries.
                {
                    // Update rule.
                    m.new_field[i * (n + 2) + j] = (1.0 / 4) * (m.current_field[(i - 1) * (n + 2) + j] + m.current_field[(i + 1) * (n + 2) + j] + m.current_field[i * (n + 2) + j - 1] + m.current_field[i * (n + 2) + j + 1]);
                }
            }
        }
        {
            SimpleTimer t("Communication Duration");
            MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
        }

        {
            SimpleTimer t("Computution Duration");
            #pragma omp parallel for
            for (int j = 1; j < m.n_loc + 1; ++j) 
            {
                m.new_field[(n + 2) + j] = (1.0 / 4) * (m.current_field[j] + m.current_field[2 * (n + 2) + j] + m.current_field[(n + 2) + j - 1] + m.current_field[(n + 2) + j + 1]);
                m.new_field[(m.n_loc) * (n + 2) + j] = (1.0 / 4) * (m.current_field[(m.n_loc - 1) * (n + 2) + j] + m.current_field[(m.n_loc + 1) * (n + 2) + j] + m.current_field[(m.n_loc) * (n + 2) + j - 1] + m.current_field[(m.n_loc) * (n + 2) + j + 1]);
            }
        }
        // Swap the current field values with the new field values before the next iteration
        std::swap(m.new_field, m.current_field);

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

/**
 * Parallel non-blokcing solver.
 */
template <typename T>
void Solver<T>::jacobi_parallel_nonblocking()
{
    int n = m.size;
    /**
     * Send last row to the next process,
     * and the first row to the previous processes.
     */
    for (int k = 0; k < max_steps; k++) // Iterate over the total number of steps
    {

        MPI_Request requests[4]; // Two sends and two receives

        // Send to the next process and receive from the previous process
        MPI_Irecv(m.current_field.data(), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == 0 ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, &requests[0]);
        MPI_Irecv(m.current_field.data() + (m.size + 2) * (m.n_loc + 1), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == m.npes - 1 ? MPI_PROC_NULL : m.rank + 1, 0, m.comm, &requests[1]);
        // Send to the previous process and receive from the next process

        MPI_Isend(m.current_field.data() + (m.size + 2) * m.n_loc, (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == m.npes - 1 ? MPI_PROC_NULL : m.rank + 1, 0, m.comm, &requests[2]);
        MPI_Isend(m.current_field.data() + (m.size + 2), (m.size + 2) * sizeof(T), MPI_BYTE,
                  m.rank == 0 ? MPI_PROC_NULL : m.rank - 1, 0, m.comm, &requests[3]);
        {
            SimpleTimer t("Computution Duration");

            for (int i = 2; i < m.n_loc; ++i) // Iterate the rows excluding the boundaries.
            {
                for (int j = 1; j < n + 1; ++j) // Iterate the columns excluding the boundaries.
                {
                    // Update rule.
                    m.new_field[i * (n + 2) + j] = (1.0 / 4) * (m.current_field[(i - 1) * (n + 2) + j] + m.current_field[(i + 1) * (n + 2) + j] + m.current_field[i * (n + 2) + j - 1] + m.current_field[i * (n + 2) + j + 1]);
                }
            }
        }
        {
            SimpleTimer t("Communication Duration");
            MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
        }
        {
            SimpleTimer t("Computution Duration");

            for (int j = 1; j < m.n_loc + 1; ++j) // Iterate the rows excluding the boundaries.
            {
                m.new_field[(n + 2) + j] = (1.0 / 4) * (m.current_field[j] + m.current_field[2 * (n + 2) + j] + m.current_field[(n + 2) + j - 1] + m.current_field[(n + 2) + j + 1]);
                m.new_field[(m.n_loc) * (n + 2) + j] = (1.0 / 4) * (m.current_field[(m.n_loc - 1) * (n + 2) + j] + m.current_field[(m.n_loc + 1) * (n + 2) + j] + m.current_field[(m.n_loc) * (n + 2) + j - 1] + m.current_field[(m.n_loc) * (n + 2) + j + 1]);
            }
        }
        // Swap the current field values with the new field values before the next iteration
        std::swap(m.new_field, m.current_field);

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
/*Jacobi hybrid*/
template <typename T>
void Solver<T>::jacobi_parallel_blocking()
{
    int n = m.size;
    int counter = 1;
    /**
     * Send last row to the next process,
     * and the first row to the previous processes.
     */
    for (int k = 0; k < max_steps; k++) // Iterate over the total number of steps
    {

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

        {
            SimpleTimer t("Computution Duration");

            for (int i = 1; i < m.n_loc + 1; ++i) // Iterate the rows excluding the boundaries.
            {
                for (int j = 1; j < n + 1; ++j) // Iterate the columns excluding the boundaries.
                {
                    // Update rule.
                    m.new_field[i * (n + 2) + j] = (1.0 / 4) * (m.current_field[(i - 1) * (n + 2) + j] + m.current_field[(i + 1) * (n + 2) + j] + m.current_field[i * (n + 2) + j - 1] + m.current_field[i * (n + 2) + j + 1]);
                }
            }
        }
        // Swap the current field values with the new field values before the next iteration
        std::swap(m.new_field, m.current_field);

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
void Solver<T>::print_results(Mesh<T> &m, int &counter)
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
#endif // JACOBI_H