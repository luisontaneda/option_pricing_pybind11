project(option_pricing_code)
cmake_minimum_required(VERSION 3.10)

add_subdirectory(thirdparty/pybind11)

pybind11_add_module(option_pricing 
                    wrappers.cpp)