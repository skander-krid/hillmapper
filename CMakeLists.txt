cmake_minimum_required(VERSION 3.15)
project(hillmapper LANGUAGES CXX)

# Require C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Fetch pybind11 automatically
include(FetchContent)
FetchContent_Declare(
  pybind11
  URL https://github.com/pybind/pybind11/archive/refs/tags/v3.0.0.tar.gz
)
FetchContent_MakeAvailable(pybind11)

# Build the Python extension
pybind11_add_module(_hillmapper ${CMAKE_CURRENT_SOURCE_DIR}/src/hillmapper.cpp)

# Put the output into the package directory
set_target_properties(_hillmapper PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/hillmapper
)

# Install the compiled extension into the wheel package (hillmapper/)
install(TARGETS _hillmapper
    LIBRARY DESTINATION hillmapper  # Linux/macOS
    ARCHIVE DESTINATION hillmapper  # Static libs, if any
    RUNTIME DESTINATION hillmapper  # Windows DLLs
)
