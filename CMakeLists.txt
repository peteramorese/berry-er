cmake_minimum_required(VERSION 3.16)

project(berryer)

find_package (Eigen3 3.3.8 REQUIRED NO_MODULE)
find_package(Threads REQUIRED)
find_package(ortools REQUIRED)
#find_package(OpenMP REQUIRED)
#find_package(BLAS REQUIRED)

option(BRYR_USE_MKL "Use Intel oneMKL linear algebra operations" OFF)
if(BRYR_USE_MKL)
    message("Using Intel oneMKL")
    set(CMAKE_PREFIX_PATH "/opt/intel/oneapi/mkl/2024.2/lib/cmake/mkl" ${CMAKE_PREFIX_PATH})
    find_package(MKL CONFIG REQUIRED)
    add_definitions(-DEIGEN_USE_MKL_ALL)
endif()

option(BRYR_USE_CUBLAS "Use cuBLAS (NVIDIA CUDA) matrix multiplication" OFF)
if(BRYR_USE_CUBLAS)
    message("Supporting NVIDIA cuBLAS")
	find_package(CUDA REQUIRED)
	find_package(CUDAToolkit REQUIRED)
    add_definitions(-DBRYR_USE_CUBLAS)
    enable_language(CUDA)

    option(CUDA_ARCHITECTURES "Specify CUDA architectures (default: 75 (GeForce RTX 2070))" "75")
    string(REPLACE " " ";" CUDA_ARCHITECTURES_LIST "${CUDA_ARCHITECTURES}")
    set(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCHITECTURES_LIST})
endif()

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
    #OpenMP::OpenMP_CXX
    #${BLAS_LIBRARIES}
)

if(BRYR_USE_MKL)
    set(BRYR_DEPENDENCY_INCLUDE_DIRS ${BRYR_DEPENDENCY_INCLUDE_DIRS} $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    set(BRYR_DEPENDENCY_LIBS ${BRYR_DEPENDENCY_LIBS} $<LINK_ONLY:MKL::MKL>)
endif()

if(BRYR_USE_CUBLAS)
    set(BRYR_DEPENDENCY_INCLUDE_DIRS ${BRYR_DEPENDENCY_INCLUDE_DIRS} ${CUDA_INCLUDE_DIRS})
    set(BRYR_DEPENDENCY_LIBS ${BRYR_DEPENDENCY_LIBS} ${CUDA_CUBLAS_LIBRARIES} CUDA::cudart)
endif()

set(BRYR_INCLUDE_DIRS 
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/include/impl
    ${BRYR_DEPENDENCY_INCLUDE_DIRS}
)

if(BRYR_USE_CUBLAS)
    add_subdirectory(cublas_src)
endif()

add_subdirectory(src)
