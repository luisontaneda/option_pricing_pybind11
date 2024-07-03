cmake_minimum_required(VERSION 3.10)

project(option_pricing_code)


# Add include directories
include_directories(
    ${CMAKE_SOURCE_DIR}/thirdparty/eigen
    ${CMAKE_SOURCE_DIR}/pricing_algorithms/COS_method
    ${CMAKE_SOURCE_DIR}/pricing_algorithms/montecarlo_method
)

add_subdirectory(thirdparty/pybind11)


pybind11_add_module(option_pricing 
                    wrappers.cpp)


target_sources(option_pricing PRIVATE
     pricing_algorithms/COS_method/cos_pricing.h
     pricing_algorithms/montecarlo_method/stoch_proccess.h
 )
