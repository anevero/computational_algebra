cmake_minimum_required(VERSION 3.14)
project(cma)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

if (MINGW)
    set(THREADS_OPTION "-mthreads")
else (MINGW)
    set(THREADS_OPTION "-pthread")
endif (MINGW)

set(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} "-no-pie ${THREADS_OPTION}")

add_executable(cma main.cpp Matrix.h
        Tests/Tests.h Tests/Tests.cpp
        Tasks/Task1.h Tasks/Task2.h Tasks/Task3.h
        Utils/Utils.cpp Utils/Utils.h Utils/DataOutput.cpp Utils/DataOutput.h
        Utils/thread_pool.cpp Utils/thread_pool.h)
