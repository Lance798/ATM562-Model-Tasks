cmake_minimum_required(VERSION 4.2)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 23)

project(ModelTasks)

set(SOURCES
    Src/main.cpp
    Src/MpiHandler.cpp
)

add_compile_options(-stdlib=libc++)
add_link_options(-stdlib=libc++)

set(CMAKE_BUILD_TYPE Debug)


find_package(MPI REQUIRED)

add_executable(${PROJECT_NAME} ${SOURCES})
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Inc)
target_link_libraries(${PROJECT_NAME} PRIVATE MPI::MPI_CXX)
target_compile_options(${PROJECT_NAME} PRIVATE -O3)
target_compile_options(${PROJECT_NAME} PRIVATE -march=native)
target_compile_options(${PROJECT_NAME} PRIVATE -Werror)
