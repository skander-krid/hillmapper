cmake_minimum_required(VERSION 3.15)
project(hillmapper LANGUAGES CXX)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)

FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v2.13.4.tar.gz
)
FetchContent_MakeAvailable(pybind11)

pybind11_add_module(_hillmapper ${CMAKE_CURRENT_SOURCE_DIR}/src/hillmapper.cpp)

find_package(pybind11 REQUIRED PATHS ${pybind11_DIR})

set_target_properties(_hillmapper PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/hillmapper
)
