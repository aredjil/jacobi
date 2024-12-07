#pragma once
#include <chrono>
#include <string>
#include <map>
#include <iostream>
#include <cmath>
#include <vector>
#include <mpi.h>
#include <cblas.h>

using units = std::chrono::microseconds;

template <typename T>
struct Timer_data
{
  T time;
  int calls;
};

template <typename T>
class CSimple_timer
{
public:
  struct Timer_data<T> timer;
  std::chrono::time_point<std::chrono::steady_clock> start;
  std::chrono::time_point<std::chrono::steady_clock> end;
  std::string what;
  std::map<std::string, struct Timer_data<T>> *table;
  bool call_destructor = true;

  CSimple_timer(const std::string &timewhat0, std::map<std::string, struct Timer_data<T>> &table)
  {
    this->start = std::chrono::steady_clock::now();
    this->what = timewhat0;
    this->table = &table;

    if (table.find(timewhat0) == table.end())
    {
      this->timer.calls = 0;
      table[timewhat0] = this->timer;
      table[timewhat0].time = std::chrono::duration_cast<T>(this->start - this->start);
    }
  }

  static void mpi_send_map(std::map<std::string, struct Timer_data<T>> &table, const int &dest, const int &tag, MPI_Comm comm)
  {
    int size = table.size();

    MPI_Send(&size, 1, MPI_INT, dest, tag, comm);

    for (auto &i : table)
    {
      int string_lenght = i.first.size();
      MPI_Send(&string_lenght, 1, MPI_INT, dest, tag, comm);

      MPI_Send(i.first.c_str(), string_lenght, MPI_CHAR, dest, tag, comm);

      MPI_Send(&i.second.calls, 1, MPI_INT, dest, tag, comm);

      long long int time_casted = i.second.time.count();
      MPI_Send(&time_casted, sizeof(T), MPI_BYTE, dest, tag, comm);
    }
  }

  static void mpi_recv_map(std::map<std::string, struct Timer_data<T>> &table, const int &source, const int &tag, MPI_Comm comm)
  {
    int size;

    MPI_Recv(&size, 1, MPI_INT, source, tag, comm, MPI_STATUS_IGNORE);

    for (int i = 0; i < size; i++)
    {
      int string_lenght;
      MPI_Recv(&string_lenght, 1, MPI_INT, source, tag, comm, MPI_STATUS_IGNORE);

      std::vector<char> what_is(string_lenght + 1);
      MPI_Recv(what_is.data(), string_lenght, MPI_CHAR, source, tag, comm, MPI_STATUS_IGNORE);
      what_is[string_lenght] = '\0';

      std::string timewhat0(what_is.data());

      int other_calls;
      MPI_Recv(&other_calls, 1, MPI_INT, source, tag, comm, MPI_STATUS_IGNORE);

      long long int other_time;
      MPI_Recv(&other_time, sizeof(T), MPI_BYTE, source, tag, comm, MPI_STATUS_IGNORE);

      if (table.find(timewhat0) == table.end())
      {
        struct Timer_data<T> timer_target;
        timer_target.calls = other_calls;
        timer_target.time = T(other_time);

        table[timewhat0] = timer_target;
      }
      else
      {
        table[timewhat0].calls += other_calls;
        table[timewhat0].time += T(other_time);
      }
    }
  }

  ~CSimple_timer()
  {
    if (call_destructor)
    {
      this->end = std::chrono::steady_clock::now();
      (*table)[what].calls += 1;
      (*table)[what].time += std::chrono::duration_cast<T>(this->end - this->start);
    }
  }

  static void print_timing_results(std::map<std::string, struct Timer_data<T>> table)
  {
    std::cout << "What                     \tCalls\tTotal time" << std::endl;
    for (const auto &i : table)
    {
      std::cout << i.first << ":\t" << i.second.calls << "\t" << i.second.time.count() << std::endl;
    }
  }

  static void collect_times(std::map<std::string, struct Timer_data<T>> &new_table, std::map<std::string, struct Timer_data<T>> &ext_table)
  {
    for (const auto &i : new_table)
    {
      ext_table[i.first] = i.second;
    }
  }
};