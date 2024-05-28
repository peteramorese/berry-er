cmake_minimum_required(VERSION 3.16)

project(berry)

find_package (Eigen3 3.3.8 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)
find_package(ortools REQUIRED)

if(NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
message("Build type: ${CMAKE_BUILD_TYPE}")

if(NOT DEFINED BRYR_BUILD_YAML)
    option(BRYR_BUILD_YAML "Build the yaml-cpp dependency" ON)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Dependencies
set(BRY_BUILD_EXECUTABLES OFF)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/modules/berry)

set(BRYR_DEPENDENCY_INCLUDE_DIRS 
    ${EIGEN3_INCLUDE_DIRS}
    ${BRY_INCLUDE_DIRS}
)

set(BRYR_DEPENDENCY_LIBS
    ortools::ortools
)

set(BRYR_INCLUDE_DIRS 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/impl
    ${BRYR_DEPENDENCY_INCLUDE_DIRS}
)

add_subdirectory(src)
