cmake_minimum_required(VERSION 3.24)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

project(
    autodiff
    LANGUAGES Fortran
    VERSION 0.1.0
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STNADARD_REQUIRED TRUE)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CTest)

add_subdirectory(src)

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
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