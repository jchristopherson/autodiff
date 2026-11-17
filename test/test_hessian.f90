program test_hessian
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_hyperdual_numbers
    use autodiff_interfaces, only: hyperdual_function
    use autodiff_hessian, only: hessian
    implicit none

    real(real64), parameter :: tol = 1.0d-10
    integer :: failures
    real(real64) :: x(2), H(2,2), g(2), f

    failures = 0

    ! Mixed polynomial, trigonometric, and exponential function
    x = [1.3_real64, 0.7_real64]
    call hessian(mixed_function, x, H, g, f)
    call check_close("mixed value", f, x(1)**2 * x(2) + sin(x(2)) + exp(x(1) * x(2)))
    call check_close("mixed g1", g(1), 2.0_real64 * x(1) * x(2) + x(2) * exp(x(1) * x(2)))
    call check_close("mixed g2", g(2), x(1)**2 + cos(x(2)) + x(1) * exp(x(1) * x(2)))
    call check_close("mixed H11", H(1,1), 2.0_real64 * x(2) + x(2)**2 * exp(x(1) * x(2)))
    call check_close("mixed H12", H(1,2), 2.0_real64 * x(1) + exp(x(1) * x(2)) &
        + x(1) * x(2) * exp(x(1) * x(2)))
    call check_close("mixed H21", H(2,1), H(1,2))
    call check_close("mixed H22", H(2,2), -sin(x(2)) + x(1)**2 * exp(x(1) * x(2)))

    ! Square root, logarithm, and inverse tangent
    x = [2.0_real64, 3.0_real64]
    call hessian(root_log_function, x, H, g)
    call check_close("root-log g1", g(1), log(x(2)) / (2.0_real64 * sqrt(x(1))) &
        + 1.0_real64 / (1.0_real64 + x(1)**2))
    call check_close("root-log g2", g(2), sqrt(x(1)) / x(2))
    call check_close("root-log H11", H(1,1), -log(x(2)) / (4.0_real64 * x(1)**1.5_real64) &
        - 2.0_real64 * x(1) / (1.0_real64 + x(1)**2)**2)
    call check_close("root-log H12", H(1,2), 1.0_real64 / (2.0_real64 * sqrt(x(1)) * x(2)))
    call check_close("root-log H22", H(2,2), -sqrt(x(1)) / x(2)**2)

    ! Division and real-valued power
    x = [1.5_real64, 2.0_real64]
    call hessian(rational_function, x, H, g, f)
    call check_close("rational value", f, x(1)**3 / x(2) + x(2)**2.5_real64)
    call check_close("rational g1", g(1), 3.0_real64 * x(1)**2 / x(2))
    call check_close("rational g2", g(2), -x(1)**3 / x(2)**2 + 2.5_real64 * x(2)**1.5_real64)
    call check_close("rational H11", H(1,1), 6.0_real64 * x(1) / x(2))
    call check_close("rational H12", H(1,2), -3.0_real64 * x(1)**2 / x(2)**2)
    call check_close("rational H22", H(2,2), 2.0_real64 * x(1)**3 / x(2)**3 &
        + 3.75_real64 * sqrt(x(2)))

    if (failures > 0) then
        write(error_unit, '(A,I0,A)') "hessian tests failed: ", failures, " check(s)"
        error stop 1
    end if
    write(output_unit, '(A)') "All hessian tests passed."

contains

    function mixed_function(u, args) result(fv)
        type(hyperdual), intent(in), dimension(:) :: u
        class(*), intent(in), optional :: args
        type(hyperdual) :: fv
        fv = u(1)**2 * u(2) + sin(u(2)) + exp(u(1) * u(2))
    end function

    function root_log_function(u, args) result(fv)
        type(hyperdual), intent(in), dimension(:) :: u
        class(*), intent(in), optional :: args
        type(hyperdual) :: fv
        fv = sqrt(u(1)) * log(u(2)) + atan(u(1))
    end function

    function rational_function(u, args) result(fv)
        type(hyperdual), intent(in), dimension(:) :: u
        class(*), intent(in), optional :: args
        type(hyperdual) :: fv
        fv = u(1)**3 / u(2) + u(2)**2.5_real64
    end function

    subroutine check_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: actual, expected
        real(real64) :: scale
        scale = max(1.0_real64, abs(expected))
        if (abs(actual - expected) > tol * scale) then
            failures = failures + 1
            write(error_unit, '(A,A,ES22.14,A,ES22.14)') "FAILED: ", label, &
                actual, " /= ", expected
        end if
    end subroutine

end program
