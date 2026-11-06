cmake_minimum_required(VERSION 3.24)

project(autodiff_install_smoke LANGUAGES Fortran)

find_package(autodiff CONFIG REQUIRED)

add_executable(autodiff_install_smoke main.f90)
target_link_libraries(autodiff_install_smoke PRIVATE autodiff::autodiff)

enable_testing()
add_test(NAME autodiff.install_smoke COMMAND autodiff_install_smoke)
