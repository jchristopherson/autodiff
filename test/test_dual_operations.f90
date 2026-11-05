program test_dual_operations
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_dual_numbers, only: dual, assignment(=), operator(+), operator(-), &
        operator(*), operator(/), operator(**), operator(==), operator(/=), operator(<), &
        operator(<=), operator(>), operator(>=), abs, sqrt, exp, log, log10, sin, cos, tan, &
        asin, acos, atan, atan2, sinh, cosh, tanh, asinh, acosh, atanh, erf, erfc, min, max, sign
    implicit none

    type(dual) :: a, b, c
    real(real64) :: scalar
    integer :: checks

    checks = 0
    a = dual(2.0_real64, 1.0_real64)
    b = dual(3.0_real64, 2.0_real64)

    c = +a
    call check_dual("unary plus", c, 2.0_real64, 1.0_real64)
    c = -a
    call check_dual("unary minus", c, -2.0_real64, -1.0_real64)

    c = a + b
    call check_dual("dual addition", c, 5.0_real64, 3.0_real64)
    c = a + 4.0_real64
    call check_dual("dual plus real", c, 6.0_real64, 1.0_real64)
    c = 4.0_real64 + a
    call check_dual("real plus dual", c, 6.0_real64, 1.0_real64)
    c = a + 4
    call check_dual("dual plus integer", c, 6.0_real64, 1.0_real64)
    c = 4 + a
    call check_dual("integer plus dual", c, 6.0_real64, 1.0_real64)

    c = a - b
    call check_dual("dual subtraction", c, -1.0_real64, -1.0_real64)
    c = a - 4.0_real64
    call check_dual("dual minus real", c, -2.0_real64, 1.0_real64)
    c = 4.0_real64 - a
    call check_dual("real minus dual", c, 2.0_real64, -1.0_real64)
    c = a - 4
    call check_dual("dual minus integer", c, -2.0_real64, 1.0_real64)
    c = 4 - a
    call check_dual("integer minus dual", c, 2.0_real64, -1.0_real64)

    c = a * b
    call check_dual("dual multiplication", c, 6.0_real64, 7.0_real64)
    c = a * 4.0_real64
    call check_dual("dual times real", c, 8.0_real64, 4.0_real64)
    c = 4.0_real64 * a
    call check_dual("real times dual", c, 8.0_real64, 4.0_real64)
    c = a * 4
    call check_dual("dual times integer", c, 8.0_real64, 4.0_real64)
    c = 4 * a
    call check_dual("integer times dual", c, 8.0_real64, 4.0_real64)

    c = a / b
    call check_dual("dual division", c, 2.0_real64 / 3.0_real64, -1.0_real64 / 9.0_real64)
    c = a / 4.0_real64
    call check_dual("dual over real", c, 0.5_real64, 0.25_real64)
    c = 4.0_real64 / a
    call check_dual("real over dual", c, 2.0_real64, -1.0_real64)
    c = a / 4
    call check_dual("dual over integer", c, 0.5_real64, 0.25_real64)
    c = 4 / a
    call check_dual("integer over dual", c, 2.0_real64, -1.0_real64)

    c = a ** b
    call check_dual("dual power dual", c, 8.0_real64, 16.0_real64 * log(2.0_real64) + 12.0_real64)
    c = a ** 3.0_real64
    call check_dual("dual power real", c, 8.0_real64, 12.0_real64)
    c = 2.0_real64 ** dual(3.0_real64, 2.0_real64)
    call check_dual("real power dual", c, 8.0_real64, 16.0_real64 * log(2.0_real64))
    c = a ** 3
    call check_dual("dual power integer", c, 8.0_real64, 12.0_real64)
    c = 2 ** dual(3.0_real64, 2.0_real64)
    call check_dual("integer power dual", c, 8.0_real64, 16.0_real64 * log(2.0_real64))
    c = dual(-2.0_real64, 1.0_real64) ** 3
    call check_dual("negative base integer power", c, -8.0_real64, 12.0_real64)
    c = dual(0.0_real64, 0.0_real64) ** 0.5_real64
    call check_dual("zero base inactive power", c, 0.0_real64, 0.0_real64)

    call check_true("dual equality", a == dual(2.0_real64, -8.0_real64))
    call check_true("dual inequality", a /= b)
    call check_true("dual less-than", a < b)
    call check_true("dual less-or-equal", a <= dual(2.0_real64, 0.0_real64))
    call check_true("dual greater-than", b > a)
    call check_true("dual greater-or-equal", b >= 3.0_real64)
    call check_true("dual equal real", a == 2.0_real64)
    call check_true("real equal dual", 2.0_real64 == a)
    call check_true("dual equal integer", a == 2)
    call check_true("integer equal dual", 2 == a)
    call check_true("dual not equal real", a /= 3.0_real64)
    call check_true("real not equal dual", 3.0_real64 /= a)
    call check_true("dual not equal integer", a /= 3)
    call check_true("integer not equal dual", 3 /= a)
    call check_true("dual less than real", a < 3.0_real64)
    call check_true("real less than dual", 1.0_real64 < a)
    call check_true("dual less than integer", a < 3)
    call check_true("integer less than dual", 1 < a)
    call check_true("dual less-or-equal real", a <= 2.0_real64)
    call check_true("real less-or-equal dual", 2.0_real64 <= a)
    call check_true("dual less-or-equal integer", a <= 2)
    call check_true("integer less-or-equal dual", 2 <= a)
    call check_true("dual greater than real", b > 2.0_real64)
    call check_true("real greater than dual", 4.0_real64 > b)
    call check_true("dual greater than integer", b > 2)
    call check_true("integer greater than dual", 4 > b)
    call check_true("dual greater-or-equal real", b >= 3.0_real64)
    call check_true("real greater-or-equal dual", 3.0_real64 >= b)
    call check_true("dual greater-or-equal integer", b >= 3)
    call check_true("integer greater-or-equal dual", 3 >= b)

    c = abs(dual(-2.0_real64, 3.0_real64))
    call check_dual("absolute value", c, 2.0_real64, -3.0_real64)
    c = abs(dual(0.0_real64, 0.0_real64))
    call check_dual("absolute value at zero", c, 0.0_real64, 0.0_real64)
    c = sqrt(dual(4.0_real64, 2.0_real64))
    call check_dual("square root", c, 2.0_real64, 0.5_real64)
    c = sqrt(dual(0.0_real64, 0.0_real64))
    call check_dual("square root at zero with inactive seed", c, 0.0_real64, 0.0_real64)
    c = exp(dual(0.0_real64, 2.0_real64))
    call check_dual("exponential", c, 1.0_real64, 2.0_real64)
    c = log(dual(2.0_real64, 1.0_real64))
    call check_dual("natural logarithm", c, log(2.0_real64), 0.5_real64)
    c = log10(dual(10.0_real64, 2.0_real64))
    call check_dual("base ten logarithm", c, 1.0_real64, 2.0_real64 / (10.0_real64 * log(10.0_real64)))

    c = sin(dual(0.0_real64, 1.0_real64))
    call check_dual("sine", c, 0.0_real64, 1.0_real64)
    c = cos(dual(0.0_real64, 1.0_real64))
    call check_dual("cosine", c, 1.0_real64, 0.0_real64)
    c = tan(dual(0.0_real64, 1.0_real64))
    call check_dual("tangent", c, 0.0_real64, 1.0_real64)
    c = asin(dual(0.5_real64, 1.0_real64))
    call check_dual("inverse sine", c, asin(0.5_real64), 2.0_real64 / sqrt(3.0_real64))
    c = acos(dual(0.5_real64, 1.0_real64))
    call check_dual("inverse cosine", c, acos(0.5_real64), -2.0_real64 / sqrt(3.0_real64))
    c = atan(dual(0.5_real64, 1.0_real64))
    call check_dual("inverse tangent", c, atan(0.5_real64), 0.8_real64)
    c = atan(dual(1.0_real64, 1.0_real64), dual(1.0_real64, 0.0_real64))
    call check_dual("two-argument inverse tangent", c, 0.25_real64 * acos(-1.0_real64), 0.5_real64)
    c = atan(dual(1.0_real64, 1.0_real64), 1.0_real64)
    call check_dual("inverse tangent with real abscissa", c, 0.25_real64 * acos(-1.0_real64), 0.5_real64)
    c = atan(1.0_real64, dual(1.0_real64, 1.0_real64))
    call check_dual("inverse tangent with real ordinate", c, 0.25_real64 * acos(-1.0_real64), -0.5_real64)
    c = atan2(dual(1.0_real64, 1.0_real64), dual(1.0_real64, 0.0_real64))
    call check_dual("atan2 dual operands", c, 0.25_real64 * acos(-1.0_real64), 0.5_real64)
    c = atan2(dual(1.0_real64, 1.0_real64), 1.0_real64)
    call check_dual("atan2 real abscissa", c, 0.25_real64 * acos(-1.0_real64), 0.5_real64)
    c = atan2(1.0_real64, dual(1.0_real64, 1.0_real64))
    call check_dual("atan2 real ordinate", c, 0.25_real64 * acos(-1.0_real64), -0.5_real64)

    c = sinh(dual(0.0_real64, 1.0_real64))
    call check_dual("hyperbolic sine", c, 0.0_real64, 1.0_real64)
    c = cosh(dual(0.0_real64, 1.0_real64))
    call check_dual("hyperbolic cosine", c, 1.0_real64, 0.0_real64)
    c = tanh(dual(0.0_real64, 1.0_real64))
    call check_dual("hyperbolic tangent", c, 0.0_real64, 1.0_real64)
    c = asinh(dual(0.0_real64, 1.0_real64))
    call check_dual("inverse hyperbolic sine", c, 0.0_real64, 1.0_real64)
    c = acosh(dual(2.0_real64, 1.0_real64))
    call check_dual("inverse hyperbolic cosine", c, acosh(2.0_real64), 1.0_real64 / sqrt(3.0_real64))
    c = acosh(dual(1.0_real64, 0.0_real64))
    call check_dual("inverse hyperbolic cosine at one with inactive seed", c, 0.0_real64, 0.0_real64)
    c = atanh(dual(0.5_real64, 1.0_real64))
    call check_dual("inverse hyperbolic tangent", c, atanh(0.5_real64), 4.0_real64 / 3.0_real64)
    c = erf(dual(0.0_real64, 1.0_real64))
    call check_dual("error function", c, 0.0_real64, 2.0_real64 / sqrt(acos(-1.0_real64)))
    c = erfc(dual(0.0_real64, 1.0_real64))
    call check_dual("complementary error function", c, 1.0_real64, -2.0_real64 / sqrt(acos(-1.0_real64)))

    c = min(dual(2.0_real64, 1.0_real64), dual(3.0_real64, 2.0_real64))
    call check_dual("minimum dual operands", c, 2.0_real64, 1.0_real64)
    c = min(dual(2.0_real64, 1.0_real64), dual(2.0_real64, 2.0_real64))
    call check_dual("minimum tie chooses first", c, 2.0_real64, 1.0_real64)
    c = min(dual(2.0_real64, 1.0_real64), 1.0_real64)
    call check_dual("minimum with real", c, 1.0_real64, 0.0_real64)
    c = min(1.0_real64, dual(2.0_real64, 1.0_real64))
    call check_dual("minimum real first", c, 1.0_real64, 0.0_real64)
    c = max(dual(2.0_real64, 1.0_real64), dual(3.0_real64, 2.0_real64))
    call check_dual("maximum dual operands", c, 3.0_real64, 2.0_real64)
    c = max(dual(2.0_real64, 1.0_real64), dual(2.0_real64, 2.0_real64))
    call check_dual("maximum tie chooses first", c, 2.0_real64, 1.0_real64)
    c = max(dual(2.0_real64, 1.0_real64), 3.0_real64)
    call check_dual("maximum with real", c, 3.0_real64, 0.0_real64)
    c = max(3.0_real64, dual(2.0_real64, 1.0_real64))
    call check_dual("maximum real first", c, 3.0_real64, 0.0_real64)
    c = sign(dual(-2.0_real64, 3.0_real64), dual(-1.0_real64, 5.0_real64))
    call check_dual("sign with dual operands", c, -2.0_real64, 3.0_real64)
    c = sign(dual(-2.0_real64, 3.0_real64), 1.0_real64)
    call check_dual("sign with real second operand", c, 2.0_real64, -3.0_real64)
    c = sign(2.0_real64, dual(-1.0_real64, 5.0_real64))
    call check_dual("sign with real first operand", c, -2.0_real64, 0.0_real64)

    a = 5.0_real64
    call check_dual("real assignment clears derivative", a, 5.0_real64, 0.0_real64)
    a = 7
    call check_dual("integer assignment clears derivative", a, 7.0_real64, 0.0_real64)
    scalar = dual(4.0_real64, 9.0_real64)
    call check_close("real assignment extracts value", scalar, 4.0_real64)

    write(output_unit, '(a,i0,a)') "Passed ", checks, " dual operation checks."

contains

    subroutine check_dual(name, actual, expected_value, expected_derivative)
        character(len=*), intent(in) :: name
        type(dual), intent(in) :: actual
        real(real64), intent(in) :: expected_value, expected_derivative

        call check_close(name // " value", actual%x, expected_value)
        call check_close(name // " derivative", actual%dx, expected_derivative)
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
