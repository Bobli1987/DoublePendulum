cmake_minimum_required(VERSION 3.3)
project(DoublePendulum)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

include_directories("C:\\Program Files\\boost\\boost_1_60_0")
include_directories("C:\\Program Files\\eigen\\eigen-eigen-07105f7124f9")

set(SOURCE_FILES main.cpp)
add_executable(DoublePendulum ${SOURCE_FILES})