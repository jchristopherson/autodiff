add_executable(test_autodiff test_autodiff.f90)
target_link_libraries(test_autodiff PRIVATE autodiff::autodiff)

add_test(NAME autodiff.unit COMMAND test_autodiff)