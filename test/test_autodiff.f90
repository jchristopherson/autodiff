program test_autodiff
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_dual_numbers, only: dual, assignment(=), operator(+), operator(-), &
        operator(*), operator(/), operator(**), operator(==), operator(/=), operator(<), &
        operator(<=), operator(>), operator(>=), sqrt, exp, sin, cos, atan2, erf
    implicit none

    type(dual) :: x, y, z
    real(real64) :: value
    integer :: checks

    checks = 0

    x = dual(2.0_real64, 1.0_real64)
    y = dual(3.0_real64, 2.0_real64)

    z = x + y
    call check_close("addition value", z%x, 5.0_real64)
    call check_close("addition derivative", z%dx, 3.0_real64)

    z = x - y
    call check_close("subtraction value", z%x, -1.0_real64)
    call check_close("subtraction derivative", z%dx, -1.0_real64)

    z = x * y
    call check_close("multiplication value", z%x, 6.0_real64)
    call check_close("multiplication derivative", z%dx, 7.0_real64)

    z = x / y
    call check_close("division value", z%x, 2.0_real64 / 3.0_real64)
    call check_close("division derivative", z%dx, -1.0_real64 / 9.0_real64)

    z = x ** 3
    call check_close("integer power value", z%x, 8.0_real64)
    call check_close("integer power derivative", z%dx, 12.0_real64)

    z = sqrt(dual(4.0_real64, 2.0_real64))
    call check_close("square root value", z%x, 2.0_real64)
    call check_close("square root derivative", z%dx, 0.5_real64)

    z = exp(dual(0.0_real64, 2.0_real64))
    call check_close("exponential value", z%x, 1.0_real64)
    call check_close("exponential derivative", z%dx, 2.0_real64)

    z = sin(dual(0.0_real64, 1.0_real64))
    call check_close("sine value", z%x, 0.0_real64)
    call check_close("sine derivative", z%dx, 1.0_real64)

    z = cos(dual(0.0_real64, 1.0_real64))
    call check_close("cosine value", z%x, 1.0_real64)
    call check_close("cosine derivative", z%dx, 0.0_real64)

    z = atan2(dual(1.0_real64, 1.0_real64), dual(1.0_real64, 0.0_real64))
    call check_close("atan2 value", z%x, 0.7853981633974483_real64)
    call check_close("atan2 derivative", z%dx, 0.5_real64)

    z = erf(dual(0.0_real64, 1.0_real64))
    call check_close("error function value", z%x, 0.0_real64)
    call check_close("error function derivative", z%dx, 1.1283791670955126_real64)

    call check_true("equality ignores derivative", x == dual(2.0_real64, -1.0_real64))
    call check_true("inequality", x /= y)
    call check_true("less-than", x < y)
    call check_true("less-than-or-equal", x <= 2.0_real64)
    call check_true("greater-than", y > x)
    call check_true("greater-than-or-equal", 3.0_real64 >= y)

    x = 5.0_real64
    call check_close("constant assignment value", x%x, 5.0_real64)
    call check_close("constant assignment clears derivative", x%dx, 0.0_real64)
    value = x
    call check_close("real assignment extracts value", value, 5.0_real64)

    write(output_unit, '(a,i0,a)') "Passed ", checks, " autodiff checks."

contains

    subroutine check_close(name, actual, expected)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: actual, expected
        real(real64), parameter :: tolerance = 1.0e-12_real64
        real(real64) :: difference

        checks = checks + 1
        difference = actual - expected
        if (abs(difference) > tolerance * (1.0_real64 + abs(expected))) then
            write(error_unit, '(a,": expected ",es24.16,", got ",es24.16)') &
                trim(name), expected, actual
            error stop 1
        end if
    end subroutine

    subroutine check_true(name, condition)
        character(len=*), intent(in) :: name
        logical, intent(in) :: condition

        checks = checks + 1
        if (.not. condition) then
            write(error_unit, '(a)') trim(name) // ": expected true"
            error stop 1
        end if
    end subroutine
end program
