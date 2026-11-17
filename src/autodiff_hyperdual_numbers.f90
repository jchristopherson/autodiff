module autodiff_hyperdual_numbers
    !! Second-order forward-mode hyper-dual numbers and overloaded operations.
    use iso_fortran_env
    implicit none
    private
    public :: hyperdual
    !! Hyper-dual number type carrying a value, two first derivatives, and a mixed second derivative.
    public :: assignment(=)
    !! Assignment conversions between hyper-dual numbers and real or integer scalars.
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    !! Arithmetic operators propagating first and mixed second derivatives.
    public :: operator(==), operator(/=), operator(<), operator(<=)
    !! Comparisons use values only and ignore derivative seeds.
    public :: operator(>), operator(>=)
    !! Remaining value-only ordering comparisons.
    public :: abs, sqrt, exp, log, sin, cos, tan, atan, sinh, cosh, tanh
    !! Overloaded real-valued intrinsic functions with second-order derivative propagation.

    type hyperdual
        !! A value, two directional derivatives, and their mixed second derivative.
        !! Seed `dx` and `dy` with one for each independent variable under study
        !! and zero otherwise.
        real(real64) :: x
            !! Function or independent-variable value.
        real(real64) :: dx
            !! First derivative with respect to the first seed.
        real(real64) :: dy
            !! First derivative with respect to the second seed.
        real(real64) :: dxy
            !! Mixed second derivative with respect to both seeds.
    end type

    interface assignment(=)
        !! Copy hyper-dual data, assign a constant (clearing derivatives), or extract a real value.
        module procedure :: hd_assign_hd
        module procedure :: hd_assign_real
        module procedure :: hd_assign_integer
        module procedure :: real_assign_hd
    end interface

    interface operator(+)
        !! Addition, including constant operands in either order.
        module procedure :: hd_add, hd_add_real, real_add_hd
        module procedure :: hd_add_integer, integer_add_hd
    end interface

    interface operator(-)
        !! Subtraction and unary negation, including constant operands.
        module procedure :: hd_subtract, hd_subtract_real, real_subtract_hd
        module procedure :: hd_subtract_integer, integer_subtract_hd, hd_negative
    end interface

    interface operator(*)
        !! Multiplication using the product rule to second order.
        module procedure :: hd_multiply, hd_multiply_real, real_multiply_hd
        module procedure :: hd_multiply_integer, integer_multiply_hd
    end interface

    interface operator(/)
        !! Division; the denominator must be nonzero.
        module procedure :: hd_divide, hd_divide_real, real_divide_hd
        module procedure :: hd_divide_integer, integer_divide_hd
    end interface

    interface operator(**)
        !! Powers with constant or hyper-dual exponents.  A real or hyper-dual base
        !! raised to a non-integer power must be positive.
        module procedure :: hd_pow_integer, hd_pow_real, hd_pow_hd
        module procedure :: real_pow_hd
    end interface

    interface operator(==)
        !! Value-only equality.
        module procedure :: hd_eq, hd_eq_real, real_eq_hd
    end interface

    interface operator(/=)
        !! Value-only inequality.
        module procedure :: hd_ne, hd_ne_real, real_ne_hd
    end interface

    interface operator(<)
        !! Value-only less-than comparison.
        module procedure :: hd_lt, hd_lt_real, real_lt_hd
    end interface

    interface operator(<=)
        !! Value-only less-than-or-equal comparison.
        module procedure :: hd_le, hd_le_real, real_le_hd
    end interface

    interface operator(>)
        !! Value-only greater-than comparison.
        module procedure :: hd_gt, hd_gt_real, real_gt_hd
    end interface

    interface operator(>=)
        !! Value-only greater-than-or-equal comparison.
        module procedure :: hd_ge, hd_ge_real, real_ge_hd
    end interface

    interface abs
        !! Absolute value; the derivative at zero is taken as that of the positive branch.
        module procedure :: hd_abs
    end interface

    interface sqrt
        !! Square root of a positive value.
        module procedure :: hd_sqrt
    end interface

    interface exp
        !! Exponential function.
        module procedure :: hd_exp
    end interface

    interface log
        !! Natural logarithm of a positive value.
        module procedure :: hd_log
    end interface

    interface sin
        !! Sine function.
        module procedure :: hd_sin
    end interface

    interface cos
        !! Cosine function.
        module procedure :: hd_cos
    end interface

    interface tan
        !! Tangent function.
        module procedure :: hd_tan
    end interface

    interface atan
        !! Inverse tangent function.
        module procedure :: hd_atan
    end interface

    interface sinh
        !! Hyperbolic sine function.
        module procedure :: hd_sinh
    end interface

    interface cosh
        !! Hyperbolic cosine function.
        module procedure :: hd_cosh
    end interface

    interface tanh
        !! Hyperbolic tangent function.
        module procedure :: hd_tanh
    end interface

contains

    ! Builds a constant hyper-dual number (all derivatives zero).
    pure function hd_const(c) result(h)
        real(real64), intent(in) :: c
        type(hyperdual) :: h
        h%x = c
        h%dx = 0.0_real64
        h%dy = 0.0_real64
        h%dxy = 0.0_real64
    end function

    ! Applies the chain rule for f(a) given f(a.x), f'(a.x), and f''(a.x).
    pure function hd_compose(a, f0, f1, f2) result(c)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: f0
            !! The function value at a%x.
        real(real64), intent(in) :: f1
            !! The first derivative of the function at a%x.
        real(real64), intent(in) :: f2
            !! The second derivative of the function at a%x.
        type(hyperdual) :: c
        c%x = f0
        c%dx = f1 * a%dx
        c%dy = f1 * a%dy
        c%dxy = f1 * a%dxy + f2 * a%dx * a%dy
    end function

    pure subroutine hd_assign_hd(lhs, rhs)
        type(hyperdual), intent(out) :: lhs
        type(hyperdual), intent(in) :: rhs
        lhs%x = rhs%x
        lhs%dx = rhs%dx
        lhs%dy = rhs%dy
        lhs%dxy = rhs%dxy
    end subroutine

    pure subroutine hd_assign_real(lhs, rhs)
        type(hyperdual), intent(out) :: lhs
        real(real64), intent(in) :: rhs
        lhs = hd_const(rhs)
    end subroutine

    pure subroutine hd_assign_integer(lhs, rhs)
        type(hyperdual), intent(out) :: lhs
        integer, intent(in) :: rhs
        lhs = hd_const(real(rhs, real64))
    end subroutine

    pure subroutine real_assign_hd(lhs, rhs)
        real(real64), intent(out) :: lhs
        type(hyperdual), intent(in) :: rhs
        lhs = rhs%x
    end subroutine

    ! Addition and subtraction.
    pure function hd_add(a, b) result(c)
        type(hyperdual), intent(in) :: a, b
        type(hyperdual) :: c
        c%x = a%x + b%x
        c%dx = a%dx + b%dx
        c%dy = a%dy + b%dy
        c%dxy = a%dxy + b%dxy
    end function

    pure function hd_add_real(a, b) result(c)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        type(hyperdual) :: c
        c = hd_add(a, hd_const(b))
    end function

    pure function real_add_hd(a, b) result(c)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_add(hd_const(a), b)
    end function

    pure function hd_add_integer(a, b) result(c)
        type(hyperdual), intent(in) :: a
        integer, intent(in) :: b
        type(hyperdual) :: c
        c = hd_add(a, hd_const(real(b, real64)))
    end function

    pure function integer_add_hd(a, b) result(c)
        integer, intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_add(hd_const(real(a, real64)), b)
    end function

    pure function hd_subtract(a, b) result(c)
        type(hyperdual), intent(in) :: a, b
        type(hyperdual) :: c
        c%x = a%x - b%x
        c%dx = a%dx - b%dx
        c%dy = a%dy - b%dy
        c%dxy = a%dxy - b%dxy
    end function

    pure function hd_subtract_real(a, b) result(c)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        type(hyperdual) :: c
        c = hd_subtract(a, hd_const(b))
    end function

    pure function real_subtract_hd(a, b) result(c)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_subtract(hd_const(a), b)
    end function

    pure function hd_subtract_integer(a, b) result(c)
        type(hyperdual), intent(in) :: a
        integer, intent(in) :: b
        type(hyperdual) :: c
        c = hd_subtract(a, hd_const(real(b, real64)))
    end function

    pure function integer_subtract_hd(a, b) result(c)
        integer, intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_subtract(hd_const(real(a, real64)), b)
    end function

    pure function hd_negative(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        c%x = -a%x
        c%dx = -a%dx
        c%dy = -a%dy
        c%dxy = -a%dxy
    end function

    ! Multiplication and division.
    pure function hd_multiply(a, b) result(c)
        type(hyperdual), intent(in) :: a, b
        type(hyperdual) :: c
        c%x = a%x * b%x
        c%dx = a%x * b%dx + a%dx * b%x
        c%dy = a%x * b%dy + a%dy * b%x
        c%dxy = a%x * b%dxy + a%dxy * b%x + a%dx * b%dy + a%dy * b%dx
    end function

    pure function hd_multiply_real(a, b) result(c)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        type(hyperdual) :: c
        c = hd_multiply(a, hd_const(b))
    end function

    pure function real_multiply_hd(a, b) result(c)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_multiply(hd_const(a), b)
    end function

    pure function hd_multiply_integer(a, b) result(c)
        type(hyperdual), intent(in) :: a
        integer, intent(in) :: b
        type(hyperdual) :: c
        c = hd_multiply(a, hd_const(real(b, real64)))
    end function

    pure function integer_multiply_hd(a, b) result(c)
        integer, intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_multiply(hd_const(real(a, real64)), b)
    end function

    pure function hd_reciprocal(b) result(c)
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_compose(b, 1.0_real64 / b%x, -1.0_real64 / b%x**2, 2.0_real64 / b%x**3)
    end function

    pure function hd_divide(a, b) result(c)
        type(hyperdual), intent(in) :: a, b
        type(hyperdual) :: c
        c = hd_multiply(a, hd_reciprocal(b))
    end function

    pure function hd_divide_real(a, b) result(c)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        type(hyperdual) :: c
        c = hd_divide(a, hd_const(b))
    end function

    pure function real_divide_hd(a, b) result(c)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_divide(hd_const(a), b)
    end function

    pure function hd_divide_integer(a, b) result(c)
        type(hyperdual), intent(in) :: a
        integer, intent(in) :: b
        type(hyperdual) :: c
        c = hd_divide(a, hd_const(real(b, real64)))
    end function

    pure function integer_divide_hd(a, b) result(c)
        integer, intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_divide(hd_const(real(a, real64)), b)
    end function

    ! Powers.
    pure function hd_pow_integer(a, n) result(c)
        type(hyperdual), intent(in) :: a
        integer, intent(in) :: n
        type(hyperdual) :: c
        real(real64) :: f1, f2
        f1 = 0.0_real64
        f2 = 0.0_real64
        if (n /= 0) f1 = real(n, real64) * a%x**(n - 1)
        if (n /= 0 .and. n /= 1) f2 = real(n * (n - 1), real64) * a%x**(n - 2)
        c = hd_compose(a, a%x**n, f1, f2)
    end function

    pure function hd_pow_real(a, p) result(c)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: p
        type(hyperdual) :: c
        real(real64) :: f1, f2
        f1 = 0.0_real64
        f2 = 0.0_real64
        if (p /= 0.0_real64) f1 = p * a%x**(p - 1.0_real64)
        if (p /= 0.0_real64 .and. p /= 1.0_real64) f2 = p * (p - 1.0_real64) * a%x**(p - 2.0_real64)
        c = hd_compose(a, a%x**p, f1, f2)
    end function

    pure function hd_pow_hd(a, b) result(c)
        type(hyperdual), intent(in) :: a, b
        type(hyperdual) :: c
        c = hd_exp(hd_multiply(b, hd_log(a)))
    end function

    pure function real_pow_hd(a, b) result(c)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        type(hyperdual) :: c
        c = hd_exp(b * log(a))
    end function

    ! Comparisons (values only).
    pure logical function hd_eq(a, b)
        type(hyperdual), intent(in) :: a, b
        hd_eq = a%x == b%x
    end function

    pure logical function hd_eq_real(a, b)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        hd_eq_real = a%x == b
    end function

    pure logical function real_eq_hd(a, b)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        real_eq_hd = a == b%x
    end function

    pure logical function hd_ne(a, b)
        type(hyperdual), intent(in) :: a, b
        hd_ne = a%x /= b%x
    end function

    pure logical function hd_ne_real(a, b)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        hd_ne_real = a%x /= b
    end function

    pure logical function real_ne_hd(a, b)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        real_ne_hd = a /= b%x
    end function

    pure logical function hd_lt(a, b)
        type(hyperdual), intent(in) :: a, b
        hd_lt = a%x < b%x
    end function

    pure logical function hd_lt_real(a, b)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        hd_lt_real = a%x < b
    end function

    pure logical function real_lt_hd(a, b)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        real_lt_hd = a < b%x
    end function

    pure logical function hd_le(a, b)
        type(hyperdual), intent(in) :: a, b
        hd_le = a%x <= b%x
    end function

    pure logical function hd_le_real(a, b)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        hd_le_real = a%x <= b
    end function

    pure logical function real_le_hd(a, b)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        real_le_hd = a <= b%x
    end function

    pure logical function hd_gt(a, b)
        type(hyperdual), intent(in) :: a, b
        hd_gt = a%x > b%x
    end function

    pure logical function hd_gt_real(a, b)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        hd_gt_real = a%x > b
    end function

    pure logical function real_gt_hd(a, b)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        real_gt_hd = a > b%x
    end function

    pure logical function hd_ge(a, b)
        type(hyperdual), intent(in) :: a, b
        hd_ge = a%x >= b%x
    end function

    pure logical function hd_ge_real(a, b)
        type(hyperdual), intent(in) :: a
        real(real64), intent(in) :: b
        hd_ge_real = a%x >= b
    end function

    pure logical function real_ge_hd(a, b)
        real(real64), intent(in) :: a
        type(hyperdual), intent(in) :: b
        real_ge_hd = a >= b%x
    end function

    ! Intrinsic functions.
    pure function hd_abs(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        real(real64) :: s
        s = 1.0_real64
        if (a%x < 0.0_real64) s = -1.0_real64
        c = hd_compose(a, abs(a%x), s, 0.0_real64)
    end function

    pure function hd_sqrt(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        real(real64) :: s
        s = sqrt(a%x)
        c = hd_compose(a, s, 0.5_real64 / s, -0.25_real64 / (a%x * s))
    end function

    pure function hd_exp(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        real(real64) :: e
        e = exp(a%x)
        c = hd_compose(a, e, e, e)
    end function

    pure function hd_log(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        c = hd_compose(a, log(a%x), 1.0_real64 / a%x, -1.0_real64 / a%x**2)
    end function

    pure function hd_sin(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        c = hd_compose(a, sin(a%x), cos(a%x), -sin(a%x))
    end function

    pure function hd_cos(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        c = hd_compose(a, cos(a%x), -sin(a%x), -cos(a%x))
    end function

    pure function hd_tan(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        real(real64) :: t
        t = tan(a%x)
        c = hd_compose(a, t, 1.0_real64 + t * t, 2.0_real64 * t * (1.0_real64 + t * t))
    end function

    pure function hd_atan(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        real(real64) :: d
        d = 1.0_real64 + a%x * a%x
        c = hd_compose(a, atan(a%x), 1.0_real64 / d, -2.0_real64 * a%x / d**2)
    end function

    pure function hd_sinh(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        c = hd_compose(a, sinh(a%x), cosh(a%x), sinh(a%x))
    end function

    pure function hd_cosh(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        c = hd_compose(a, cosh(a%x), sinh(a%x), cosh(a%x))
    end function

    pure function hd_tanh(a) result(c)
        type(hyperdual), intent(in) :: a
        type(hyperdual) :: c
        real(real64) :: t
        t = tanh(a%x)
        c = hd_compose(a, t, 1.0_real64 - t * t, -2.0_real64 * t * (1.0_real64 - t * t))
    end function

end module
