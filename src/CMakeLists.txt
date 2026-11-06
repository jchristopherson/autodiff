add_library(autodiff
    autodiff.f90
    autodiff_dual_numbers.f90
    autodiff_interfaces.f90
    autodiff_workspace.f90
    autodiff_gradient.f90
    autodiff_jacobian.f90
    autodiff_jvp.f90
)

add_library(autodiff::autodiff ALIAS autodiff)

option(AUTODIFF_ENABLE_IPO "Enable interprocedural optimization (LTO) when supported" ON)
if(AUTODIFF_ENABLE_IPO)
    include(CheckIPOSupported)
    check_ipo_supported(RESULT autodiff_ipo_supported OUTPUT autodiff_ipo_output LANGUAGES Fortran)
    if(autodiff_ipo_supported)
        if(CMAKE_Fortran_COMPILER_ID MATCHES "^Intel(LLVM)?$" AND NOT BUILD_SHARED_LIBS)
            message(STATUS "Disabling IPO for Intel static library builds due to linker incompatibility.")
        else()
            set_target_properties(autodiff PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
        endif()
    else()
        message(STATUS "Interprocedural optimization not supported: ${autodiff_ipo_output}")
    endif()
endif()

if(AUTODIFF_ENABLE_OPENMP)
    if(OpenMP_Fortran_FOUND)
        target_link_libraries(autodiff PUBLIC OpenMP::OpenMP_Fortran)
        target_link_options(autodiff PUBLIC "${OpenMP_Fortran_FLAGS}")
    else()
        message(STATUS "OpenMP Fortran support not found; Jacobian columns will be evaluated serially.")
    endif()
endif()

set(autodiff_module_dir "${CMAKE_CURRENT_BINARY_DIR}/modules")
if(CMAKE_CONFIGURATION_TYPES)
    set(autodiff_module_build_include_dir "${autodiff_module_dir}/$<CONFIG>")
    set(autodiff_module_install_source_dir "${autodiff_module_dir}/$<CONFIG>")
else()
    set(autodiff_module_build_include_dir "${autodiff_module_dir}")
    set(autodiff_module_install_source_dir "${autodiff_module_dir}")
endif()

set_target_properties(autodiff PROPERTIES
    Fortran_MODULE_DIRECTORY "${autodiff_module_dir}"
)

target_include_directories(autodiff PUBLIC
    "$<BUILD_INTERFACE:${autodiff_module_build_include_dir}>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/autodiff>"
)

install(
    TARGETS autodiff
    EXPORT autodiffTargets
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/autodiff"
)

install(
    DIRECTORY "${autodiff_module_install_source_dir}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/autodiff"
    FILES_MATCHING
        PATTERN "*.mod"
        PATTERN "*.smod"
)