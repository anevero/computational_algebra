cmake_minimum_required(VERSION 3.10)
project(cma)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-no-pie -pthread")

add_executable(cma main.cpp Matrix.h Polynomial.h
        Tests/Tests.h Tests/Tests.cpp
        Tasks/Task1.h Tasks/Task2.h Tasks/Task3.h Tasks/Task4.h Tasks/Task5.h
        Tasks/Task6.h Tasks/Task7.h Tasks/Task8.h
        Utils/Utils.cpp Utils/Utils.h Utils/DataOutput.cpp Utils/DataOutput.h
        Utils/ThreadPool.cpp Utils/ThreadPool.h)
