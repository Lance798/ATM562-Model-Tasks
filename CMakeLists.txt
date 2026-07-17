cmake_minimum_required(VERSION 4.2)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)

project(ModelTasks)
add_executable(ModelTasks main.cpp
    grads_io.hpp)