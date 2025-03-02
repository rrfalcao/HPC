# Minimum CMake version required
cmake_minimum_required(VERSION 3.10)

# Project name and version
project(FibonacciProject VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable target
add_executable(fibo main.cpp fibo.cpp)
