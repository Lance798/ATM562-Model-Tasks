#ifndef __MPI_HPP__
#define __MPI_HPP__

#include "main.hpp"
#include <memory>
#include <mpi.h>

constexpr int COLLECTOR = 0;

class MpiHandler {
private:
  std::vector<int> nums_arr_;
  std::vector<int> displs_;

public:
  MPI_Comm comm;
  int world_size, right, left, pos;

  int init();
  // gather all process's data
  void exchange(Grid &g);
  int gather();
  // return nx displacement, original coordinate = i + get_disp()
  int get_x_offset() { return displs_[pos] / MODEL_NZ; }

  std::unique_ptr<Grid> collect();
};

extern MpiHandler g_mpi_handler;

#endif