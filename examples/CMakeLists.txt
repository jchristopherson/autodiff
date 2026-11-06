add_executable(autodiff_jacobian_example autodiff_jacobian_example.f90)
target_link_libraries(autodiff_jacobian_example autodiff)

add_executable(autodiff_gradient_example autodiff_gradient_example.f90)
target_link_libraries(autodiff_gradient_example autodiff)

add_executable(autodiff_jvp_example autodiff_jvp_example.f90)
target_link_libraries(autodiff_jvp_example autodiff)

if (BUILD_SHARED_LIBS AND WIN32)
    add_custom_command(
        TARGET autodiff_jacobian_example
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:autodiff_jacobian_example> 
        $<TARGET_FILE_DIR:autodiff_jacobian_example>
        COMMAND_EXPAND_LISTS
    )
endif()