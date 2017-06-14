cmake_minimum_required(VERSION 3.6)
project(integral)

find_package(Threads)
set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp)
add_executable(integral ${SOURCE_FILES})
target_link_libraries(integral Threads::Threads)