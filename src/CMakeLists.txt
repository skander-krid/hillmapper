cmake_minimum_required(VERSION 3.15)
project(hillmapper LANGUAGES CXX)

find_package(pybind11 REQUIRED)

# Build the extension
pybind11_add_module(_hillmapper src/hillmapper.cpp)

# Set output dir so Python can find it
set_target_properties(_hillmapper PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/hillmapper
)