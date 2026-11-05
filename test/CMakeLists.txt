add_executable(test_autodiff test_autodiff.f90)
target_link_libraries(test_autodiff PRIVATE autodiff::autodiff)

add_test(NAME autodiff.unit COMMAND test_autodiff)

add_executable(test_dual_operations test_dual_operations.f90)
target_link_libraries(test_dual_operations PRIVATE autodiff::autodiff)

add_test(NAME autodiff.dual_operations COMMAND test_dual_operations)

add_executable(test_differentiation test_differentiation.f90)
target_link_libraries(test_differentiation PRIVATE autodiff::autodiff)

add_test(NAME autodiff.differentiation COMMAND test_differentiation)