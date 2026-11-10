cmake_minimum_required(VERSION 3.24)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

project(
    autodiff
    LANGUAGES Fortran
    VERSION 1.0.1
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STNADARD_REQUIRED TRUE)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

option(AUTODIFF_ENABLE_OPENMP "Evaluate Jacobian columns in parallel with OpenMP when available" ON)
set(AUTODIFF_CONFIG_NEEDS_OPENMP FALSE)
if(AUTODIFF_ENABLE_OPENMP)
    find_package(OpenMP COMPONENTS Fortran)
    if(OpenMP_Fortran_FOUND)
        set(AUTODIFF_CONFIG_NEEDS_OPENMP TRUE)
    endif()
endif()

add_subdirectory(src)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

option(BUILD_AUTODIFF_EXAMPLES "Build AUTODIFF examples.")
if (BUILD_AUTODIFF_EXAMPLES)
    add_subdirectory(examples)
endif()

set(autodiff_cmake_install_dir "${CMAKE_INSTALL_LIBDIR}/cmake/autodiff")

configure_package_config_file(
    cmake/autodiffConfig.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/autodiffConfig.cmake"
    INSTALL_DESTINATION "${autodiff_cmake_install_dir}"
)

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/autodiffConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY SameMajorVersion
)

install(
    EXPORT autodiffTargets
    FILE autodiffTargets.cmake
    NAMESPACE autodiff::
    DESTINATION "${autodiff_cmake_install_dir}"
)

install(
    FILES
        "${CMAKE_CURRENT_BINARY_DIR}/autodiffConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/autodiffConfigVersion.cmake"
    DESTINATION "${autodiff_cmake_install_dir}"
)