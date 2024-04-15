file(GLOB BRY_EXECUTABLES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

message("Exec files: ${BRY_EXECUTABLES}")
foreach(EXEC_FILE ${BRY_EXECUTABLES})
    get_filename_component(EXEC_NAME ${EXEC_FILE} NAME_WE)
    add_executable(${EXEC_NAME} ${EXEC_FILE})
    message("bry lib dir: ${BRY_LIBRARY_DIR}")
    target_include_directories(${EXEC_NAME} PRIVATE
        ${BRYR_INCLUDE_DIRS} 
        ${BRYR_DEPENDENCY_INCLUDE_DIRS}
    )
endforeach()
