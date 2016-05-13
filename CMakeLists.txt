cmake_minimum_required(VERSION 2.8)
project(team_notebook_Innopolis_U)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DLOCAL -O2 -Wall")

add_executable(karatsuba karatsuba.cpp)

add_executable(suffix_array suffix_array.cpp)

add_executable(cartesian cartesian.cpp)

add_executable(cut_point cut_point.cpp)

add_executable(bridges bridges.cpp)

add_executable(maxFlow maxFlow.cpp)

add_executable(minCostFlow minCostFlow.cpp)

add_executable(suffix_automaton_map suffix_automaton_map.cpp)

add_executable(suffix_automaton_arr suffix_automaton_arr.cpp)

add_executable(fft fft.cpp)

add_executable(hungarian hungarian.cpp)
