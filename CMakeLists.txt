include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG latest stable
)
FetchContent_MakeAvailable(pybind11)
cmake_minimum_required(VERSION 3.15)
project(hillmapper LANGUAGES CXX)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(pybind11 REQUIRED)

# Build the extension
pybind11_add_module(_hillmapper src/hillmapper.cpp)

# Set output dir so Python can find it
set_target_properties(_hillmapper PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/hillmapper
)
