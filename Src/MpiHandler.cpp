#include "MpiHandler.hpp"
#include "main.hpp"
#include <cstdlib>
#include <iostream>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <ostream>
#include <vector>

MpiHandler g_mpi_handler;

int MpiHandler::init() {
  MPI_Init(NULL, NULL);

  int world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_size > MODEL_NX) {
    std::cerr << "Too many process" << std::endl;
    exit(-1);
  }

  if (world_size < COLLECTOR) {
    std::cerr << "cannot assign collector" << std::endl;
    exit(-1);
  }

  int dim[1] = {0};
  MPI_Dims_create(world_size, 1, dim);
  int period[1] = {true};
  MPI_Cart_create(MPI_COMM_WORLD, 1, dim, period, true, &comm);
  if (comm == MPI_COMM_NULL) {
    std::cerr << "Cannot create communicator." << std::endl;
    exit(-1);
  }
  nums_arr_.resize(world_size);
  displs_.resize(world_size);

  for (int i = 0; i < world_size; i++) {
    // don't send lateral boundary
    nums_arr_[i] = MODEL_NZ * (MODEL_NX / g_mpi_handler.world_size +
                               ((MODEL_NX % g_mpi_handler.world_size) > i));
    displs_[i] = (i == 0) ? 0 : displs_[i - 1] + nums_arr_[i - 1];
  }

  MPI_Comm_rank(comm, &world_rank);
  MPI_Cart_coords(comm, world_rank, 1, &pos);
  MPI_Cart_shift(comm, 0, 1, &left, &right);

  if (left == MPI_PROC_NULL || right == MPI_PROC_NULL) {
    std::cerr << "Period process boundary doesn't work!" << std::endl;
    exit(-1);
  }
  return pos;
}

void MpiHandler::exchange(Grid &g) {
  // send grid[1] to left and grid[nx-2] to right
  // receive from left save to grid[0] and
  // receive from right save to grid[nx-1]
  // tag is according to VarId
  int tag = 0; // VarId

  for (auto &var : g.vars) {
    MPI_Request requests[4];
    MPI_Isend(var[1].data(), MODEL_NZ, MPI_DOUBLE, left, tag << 10 | 1, comm,
              &requests[0]);
    MPI_Isend(var[grid_nx - 2].data(), MODEL_NZ, MPI_DOUBLE, right,
              tag << 10 | 2, comm, &requests[1]);
    MPI_Irecv(var[0].data(), MODEL_NZ, MPI_DOUBLE, left, tag << 10 | 2, comm,
              &requests[2]);

    MPI_Irecv(var[grid_nx - 1].data(), MODEL_NZ, MPI_DOUBLE, right,
              tag << 10 | 1, comm, &requests[3]);

    MPI_Waitall(4, requests, MPI_STATUSES_IGNORE);
    tag++;
  }
}

std::unique_ptr<Grid> MpiHandler::collect() {
  // calculate other's nx
  std::unique_ptr<Grid> all_data = nullptr;
  if (pos == COLLECTOR) {
    // receive from others or calcuate by myself
    if (std::reduce(nums_arr_.begin(), nums_arr_.end()) !=
        MODEL_NX * MODEL_NZ) {
      std::cerr << std::format("sum of all_nx is not correct, sum({}) != {}",
                               nums_arr_, MODEL_NX)
                << std::endl;
      exit(-1);
    }
    all_data = std::make_unique<Grid>();
    all_data->resize_nx(MODEL_NX);
  }

  int id = 0;
  for (auto &var : now->vars) {
    if (pos != COLLECTOR) {
      MPI_Gatherv(var[1].data(), (grid_nx - 2) * MODEL_NZ, MPI_DOUBLE, nullptr,
                  nullptr, nullptr, MPI_DOUBLE, COLLECTOR, comm);
    } else {
      MPI_Gatherv(var[1].data(), (grid_nx - 2) * MODEL_NZ, MPI_DOUBLE,
                  all_data->vars[id].data(), nums_arr_.data(), displs_.data(),
                  MPI_DOUBLE, COLLECTOR, comm);
    }
    id++;
  }
  return all_data;
}
