cmake_minimum_required(VERSION 3.10)
project(cma)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mavx2")
set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-no-pie -pthread")

add_executable(cma main.cpp matrix.h polynomial.h
        tests/tests.cpp
        utils/data_output.cpp
        utils/utils.cpp
        utils/thread_pool.cpp)
