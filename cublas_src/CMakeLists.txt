# Calibrate
set(SZ_THR_FILEPATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/size_threshold")
if(EXISTS ${SZ_THR_FILEPATH})
    message(STATUS "Calibration file found: ${SZ_THR_FILEPATH}")

    file(READ ${SZ_THR_FILEPATH} SIZE_THRESHOLD)
    message("   [Size threshold: ${SIZE_THRESHOLD}]")
else()
    message(WARNING "No size threshold file found. Consider calibrating (./calibrateMatMul) and re-running CMake for best performance")
    set(SIZE_THRESHOLD "1600")
endif()

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/SizeThreshold.h.in
    ${CMAKE_BINARY_DIR}/include/SizeThreshold.h
    @ONLY
)

# Build
file(GLOB BRYR_CUBLAS_LIB RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cu")

add_library(bryrCUBLAS STATIC ${BRYR_CUBLAS_LIB})

set_target_properties(bryrCUBLAS PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_STANDARD 14
)

target_include_directories(bryrCUBLAS PRIVATE ${BRYR_INCLUDE_DIRS})
target_link_libraries(bryrCUBLAS ${BRYR_DEPENDENCY_LIBS})
target_include_directories(bryrCUBLAS PRIVATE ${CMAKE_BINARY_DIR}/include)

