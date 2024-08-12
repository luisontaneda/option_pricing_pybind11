cmake_minimum_required(VERSION 3.10)

project(option_pricing_code)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Add include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/thirdparty/eigen
)

# Add the pybind11 subdirectory
add_subdirectory(external/pybind11)

# Source files for the module
set(SOURCE_FILES
    wrappers.cpp
    pricing_algorithms/COS_method/European/stoch_process.cpp
    pricing_algorithms/montecarlo_method/European/stoch_process.cpp
    pricing_algorithms/COS_method/American/stoch_process.cpp
    pricing_algorithms/montecarlo_method/American/stoch_process.cpp
)

# Create the pybind11 module
pybind11_add_module(option_pricing ${SOURCE_FILES})

# Include directories for the module
target_include_directories(option_pricing PRIVATE
    ${CMAKE_SOURCE_DIR}/pricing_algorithms/COS_method
    ${CMAKE_SOURCE_DIR}/pricing_algorithms/montecarlo_method
)

# Link libraries (if any additional linking is needed, such as with Boost)
# target_link_libraries(option_pricing PRIVATE ${YOUR_LIBRARIES_HERE})