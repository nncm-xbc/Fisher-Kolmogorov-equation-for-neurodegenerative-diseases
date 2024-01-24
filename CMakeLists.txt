cmake_minimum_required(VERSION 3.12.0)
project(FKneuro LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(fish-exe src/main.cpp src/HeatNonLinear.cpp)
deal_ii_setup_target(fish-exe)