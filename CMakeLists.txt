cmake_minimum_required(VERSION 3.10)

# Set the project name
project(obj2cut)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(obj2cut obj2cut.cc)

# Include directories
include_directories(include)

# Link libraries
target_link_libraries(obj2cut ${CMAKE_THREAD_LIBS_INIT})