file(GLOB BRY_EXECUTABLES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

foreach(EXEC_FILE ${BRY_EXECUTABLES})
    get_filename_component(EXEC_NAME ${EXEC_FILE} NAME_WE)
    add_executable(${EXEC_NAME} ${EXEC_FILE})
    target_include_directories(${EXEC_NAME} PRIVATE
        ${BRYR_INCLUDE_DIRS} 
    )
    target_link_libraries(${EXEC_NAME} ${BRYR_DEPENDENCY_LIBS})
endforeach()
