program test_hyperdual_operations
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_hyperdual_numbers, only: hyperdual, abs, sqrt, exp, log, log10, sin, cos, tan, &
        asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, erf, erfc, min, max, sign
    implicit none

    type(hyperdual) :: a, b, c
    real(real64), parameter :: pi = acos(-1.0_real64)
    real(real64) :: d
    integer :: checks

    checks = 0
    a = hyperdual(0.5_real64, 2.0_real64, -1.0_real64, 3.0_real64)

    c = abs(hyperdual(-0.5_real64, 2.0_real64, -1.0_real64, 3.0_real64))
    call check_composition("absolute value", c, 0.5_real64, -1.0_real64, 0.0_real64)

    c = sqrt(hyperdual(4.0_real64, 2.0_real64, -1.0_real64, 3.0_real64))
    call check_composition("square root", c, 2.0_real64, 0.25_real64, -0.03125_real64)
    c = exp(a)
    call check_composition("exponential", c, exp(0.5_real64), exp(0.5_real64), exp(0.5_real64))
    c = log(a)
    call check_composition("natural logarithm", c, log(0.5_real64), 2.0_real64, -4.0_real64)
    c = log10(a)
    d = log(10.0_real64)
    call check_composition("base-ten logarithm", c, log10(0.5_real64), &
        1.0_real64 / (0.5_real64 * d), -1.0_real64 / (0.25_real64 * d))

    c = sin(a)
    call check_composition("sine", c, sin(0.5_real64), cos(0.5_real64), -sin(0.5_real64))
    c = cos(a)
    call check_composition("cosine", c, cos(0.5_real64), -sin(0.5_real64), -cos(0.5_real64))
    c = tan(a)
    d = tan(0.5_real64)
    call check_composition("tangent", c, d, 1.0_real64 + d * d, 2.0_real64 * d * (1.0_real64 + d * d))

    c = asin(a)
    d = 1.0_real64 - 0.5_real64**2
    call check_composition("inverse sine", c, asin(0.5_real64), 1.0_real64 / sqrt(d), &
        0.5_real64 / d**1.5_real64)
    c = acos(a)
    call check_composition("inverse cosine", c, acos(0.5_real64), -1.0_real64 / sqrt(d), &
        -0.5_real64 / d**1.5_real64)
    c = atan(a)
    d = 1.0_real64 + 0.5_real64**2
    call check_composition("inverse tangent", c, atan(0.5_real64), 1.0_real64 / d, &
        -1.0_real64 / d**2)

    b = hyperdual(2.0_real64, 0.0_real64, 1.0_real64, 0.0_real64)
    a = hyperdual(1.0_real64, 1.0_real64, 0.0_real64, 0.0_real64)
    c = atan(a, b)
    call check_hyperdual("two-argument inverse tangent", c, atan2(1.0_real64, 2.0_real64), &
        0.4_real64, -0.2_real64, -3.0_real64 / 25.0_real64)
    c = atan(a, 2.0_real64)
    call check_hyperdual("inverse tangent with real abscissa", c, atan2(1.0_real64, 2.0_real64), &
        0.4_real64, 0.0_real64, 0.0_real64)
    c = atan(1.0_real64, b)
    call check_hyperdual("inverse tangent with real ordinate", c, atan2(1.0_real64, 2.0_real64), &
        0.0_real64, -0.2_real64, 0.0_real64)
    c = atan2(a, b)
    call check_hyperdual("atan2 hyperdual operands", c, atan2(1.0_real64, 2.0_real64), &
        0.4_real64, -0.2_real64, -3.0_real64 / 25.0_real64)
    c = atan2(a, 2.0_real64)
    call check_hyperdual("atan2 real abscissa", c, atan2(1.0_real64, 2.0_real64), &
        0.4_real64, 0.0_real64, 0.0_real64)
    c = atan2(1.0_real64, b)
    call check_hyperdual("atan2 real ordinate", c, atan2(1.0_real64, 2.0_real64), &
        0.0_real64, -0.2_real64, 0.0_real64)

    a = hyperdual(0.5_real64, 2.0_real64, -1.0_real64, 3.0_real64)
    c = sinh(a)
    call check_composition("hyperbolic sine", c, sinh(0.5_real64), cosh(0.5_real64), sinh(0.5_real64))
    c = cosh(a)
    call check_composition("hyperbolic cosine", c, cosh(0.5_real64), sinh(0.5_real64), cosh(0.5_real64))
    c = tanh(a)
    d = tanh(0.5_real64)
    call check_composition("hyperbolic tangent", c, d, 1.0_real64 - d * d, &
        -2.0_real64 * d * (1.0_real64 - d * d))
    c = asinh(a)
    d = 1.0_real64 + 0.5_real64**2
    call check_composition("inverse hyperbolic sine", c, asinh(0.5_real64), 1.0_real64 / sqrt(d), &
        -0.5_real64 / d**1.5_real64)
    c = acosh(hyperdual(2.0_real64, 2.0_real64, -1.0_real64, 3.0_real64))
    d = 2.0_real64**2 - 1.0_real64
    call check_composition("inverse hyperbolic cosine", c, acosh(2.0_real64), 1.0_real64 / sqrt(d), &
        -2.0_real64 / d**1.5_real64)
    c = atanh(a)
    d = 1.0_real64 - 0.5_real64**2
    call check_composition("inverse hyperbolic tangent", c, atanh(0.5_real64), 1.0_real64 / d, &
        2.0_real64 * 0.5_real64 / d**2)
    c = erf(a)
    d = 2.0_real64 / sqrt(pi) * exp(-0.5_real64**2)
    call check_composition("error function", c, erf(0.5_real64), d, -2.0_real64 * 0.5_real64 * d)
    c = erfc(a)
    d = -2.0_real64 / sqrt(pi) * exp(-0.5_real64**2)
    call check_composition("complementary error function", c, erfc(0.5_real64), &
        d, -2.0_real64 * 0.5_real64 * d)

    a = hyperdual(2.0_real64, 1.0_real64, 2.0_real64, 3.0_real64)
    b = hyperdual(3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64)
    c = min(a, b)
    call check_hyperdual("minimum hyperdual operands", c, 2.0_real64, 1.0_real64, 2.0_real64, 3.0_real64)
    c = min(a, 1.0_real64)
    call check_hyperdual("minimum with real", c, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64)
    c = min(1.0_real64, a)
    call check_hyperdual("minimum real first", c, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64)
    c = min(a, hyperdual(2.0_real64, 4.0_real64, 5.0_real64, 6.0_real64))
    call check_hyperdual("minimum tie chooses first", c, 2.0_real64, 1.0_real64, 2.0_real64, 3.0_real64)
    c = max(a, b)
    call check_hyperdual("maximum hyperdual operands", c, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64)
    c = max(a, 3.0_real64)
    call check_hyperdual("maximum with real", c, 3.0_real64, 0.0_real64, 0.0_real64, 0.0_real64)
    c = max(3.0_real64, a)
    call check_hyperdual("maximum real first", c, 3.0_real64, 0.0_real64, 0.0_real64, 0.0_real64)
    c = max(a, hyperdual(2.0_real64, 4.0_real64, 5.0_real64, 6.0_real64))
    call check_hyperdual("maximum tie chooses first", c, 2.0_real64, 1.0_real64, 2.0_real64, 3.0_real64)

    a = hyperdual(-2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64)
    b = hyperdual(-1.0_real64, 7.0_real64, 8.0_real64, 9.0_real64)
    c = sign(a, b)
    call check_hyperdual("sign with hyperdual operands", c, -2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64)
    c = sign(a, 1.0_real64)
    call check_hyperdual("sign with real second operand", c, 2.0_real64, -3.0_real64, -4.0_real64, -5.0_real64)
    c = sign(2.0_real64, b)
    call check_hyperdual("sign with real first operand", c, -2.0_real64, 0.0_real64, 0.0_real64, 0.0_real64)

    write(output_unit, '(a,i0,a)') "Passed ", checks, " hyperdual intrinsic checks."

contains

    subroutine check_composition(name, actual, expected_value, first, second)
        character(len=*), intent(in) :: name
        type(hyperdual), intent(in) :: actual
        real(real64), intent(in) :: expected_value, first, second

        call check_hyperdual(name, actual, expected_value, 2.0_real64 * first, &
            -first, 3.0_real64 * first - 2.0_real64 * second)
    end subroutine

    subroutine check_hyperdual(name, actual, expected_value, expected_dx, expected_dy, expected_dxy)
        character(len=*), intent(in) :: name
        type(hyperdual), intent(in) :: actual
        real(real64), intent(in) :: expected_value, expected_dx, expected_dy, expected_dxy

        call check_close(name // " value", actual%x, expected_value)
        call check_close(name // " dx", actual%dx, expected_dx)
        call check_close(name // " dy", actual%dy, expected_dy)
        call check_close(name // " dxy", actual%dxy, expected_dxy)
    end subroutine

    subroutine check_close(name, actual, expected)
        character(len=*), intent(in) :: name
        real(real64), intent(in) :: actual, expected
        real(real64), parameter :: tolerance = 2.0e-12_real64

        checks = checks + 1
        if (abs(actual - expected) > tolerance * (1.0_real64 + abs(expected))) then
            write(error_unit, '(a,": expected ",es24.16,", got ",es24.16)') &
                trim(name), expected, actual
            error stop 1
        end if
    end subroutine

end program
