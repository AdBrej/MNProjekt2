cmake_minimum_required(VERSION 3.7)
project(MN2)

set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILES main.cpp constValue.h tmp.cpp)
add_executable(MN2 ${SOURCE_FILES})