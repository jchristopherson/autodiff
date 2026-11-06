module autodiff_dual_numbers
    !! Forward-mode dual numbers and overloaded arithmetic and intrinsic operations.
    use iso_fortran_env
    !! Supplies the real64 kind used for values and directional derivatives.
    implicit none
    private
    public :: dual
    !! Dual-number value type; each instance stores a value and one directional derivative.
    public :: assignment(=)
    !! Assignment conversions between dual numbers and real or integer scalars.
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    !! Arithmetic operators propagate derivatives using the corresponding calculus rules.
    public :: operator(==), operator(/=), operator(<), operator(<=)
    !! Equality and ordering comparisons use values only and ignore derivative seeds.
    public :: operator(>), operator(>=)
    !! Remaining value-only ordering comparisons.
    public :: abs, sqrt, exp, log, log10, sin, cos, tan, asin, acos, atan
    !! Overloaded real-valued intrinsic functions with derivative propagation.
    public :: atan2, sinh, cosh, tanh, asinh, acosh, atanh, erf, erfc
    !! Additional inverse trigonometric, hyperbolic, and error functions.
    public :: min, max, sign
    !! Selection and sign intrinsics propagate the derivative of the selected or magnitude operand.

    type dual
        !! A value and its directional derivative. Seed dx with one for an
        !! independent variable and zero for constants.
        real(real64) :: x
            !! Function or independent-variable value.
        real(real64) :: dx
            !! Directional derivative with respect to the current seed.
    end type

    interface assignment(=)
        !! Copy AD data, assign a constant (clearing dx), or extract a real value.
        module procedure :: ad_assign_1
        module procedure :: ad_assign_2
        module procedure :: ad_assign_3
        module procedure :: ad_assign_integer
    end interface

    interface operator(+)
        !! Unary identity and addition, including constant operands in either order.
        module procedure :: ad_add, ad_add_real, duald_ad
        module procedure :: ad_add_integer, integer_add_ad, ad_positive
    end interface
    interface operator(-)
        !! Unary negation and subtraction, including constant operands.
        module procedure :: ad_subtract, ad_subtract_real, real_subtract_ad
        module procedure :: ad_subtract_integer, integer_subtract_ad, ad_negative
    end interface
    interface operator(*)
        !! Multiplication using the product rule.
        module procedure :: ad_multiply, ad_multiply_real, real_multiply_ad
        module procedure :: ad_multiply_integer, integer_multiply_ad
    end interface
    interface operator(/)
        !! Division using the quotient rule; the denominator must be nonzero.
        module procedure :: ad_divide, ad_divide_real, real_divide_ad
        module procedure :: ad_divide_integer, integer_divide_ad
    end interface
    interface operator(**)
        !! Exponentiation with constant or variable base/exponent.
        !! Variable exponents require a positive base, except that a zero constant
        !! base with a positive exponent has zero derivative. Constant integral
        !! exponents also support negative bases.
        module procedure :: ad_power, ad_power_real, real_power_ad
        module procedure :: ad_power_integer, integer_power_ad
    end interface
    interface operator(==)
        !! Equality of values, independent of derivative seeds.
        module procedure :: ad_equal, ad_equal_real, real_equal_ad
        module procedure :: ad_equal_integer, integer_equal_ad
    end interface
    interface operator(/=)
        !! Inequality of values, independent of derivative seeds.
        module procedure :: ad_not_equal, ad_not_equal_real, real_not_equal_ad
        module procedure :: ad_not_equal_integer, integer_not_equal_ad
    end interface
    interface operator(<)
        !! Strict less-than comparison of values.
        module procedure :: ad_less, ad_less_real, real_less_ad
        module procedure :: ad_less_integer, integer_less_ad
    end interface
    interface operator(<=)
        !! Less-than-or-equal comparison of values.
        module procedure :: ad_less_equal, ad_less_equal_real, real_less_equal_ad
        module procedure :: ad_less_equal_integer, integer_less_equal_ad
    end interface
    interface operator(>)
        !! Strict greater-than comparison of values.
        module procedure :: ad_greater, ad_greater_real, real_greater_ad
        module procedure :: ad_greater_integer, integer_greater_ad
    end interface
    interface operator(>=)
        !! Greater-than-or-equal comparison of values.
        module procedure :: ad_greater_equal, ad_greater_equal_real, real_greater_equal_ad
        module procedure :: ad_greater_equal_integer, integer_greater_equal_ad
    end interface
    interface abs
        !! Absolute value; the derivative at zero is defined to be zero.
        module procedure :: ad_abs
    end interface
    interface sqrt
        !! Square root for nonnegative values; its active derivative is singular at zero.
        module procedure :: ad_sqrt
    end interface
    interface exp
        !! Natural exponential.
        module procedure :: ad_exp
    end interface
    interface log
        !! Natural logarithm for positive values.
        module procedure :: ad_log
    end interface
    interface log10
        !! Base-ten logarithm for positive values.
        module procedure :: ad_log10
    end interface
    interface sin
        !! Sine, with arguments in radians.
        module procedure :: ad_sin
    end interface
    interface cos
        !! Cosine, with arguments in radians.
        module procedure :: ad_cos
    end interface
    interface tan
        !! Tangent, with arguments in radians and away from its poles.
        module procedure :: ad_tan
    end interface
    interface asin
        !! Inverse sine in radians; the active derivative is singular at +/-1.
        module procedure :: ad_asin
    end interface
    interface acos
        !! Inverse cosine in radians; the active derivative is singular at +/-1.
        module procedure :: ad_acos
    end interface
    interface atan
        !! One- or two-argument inverse tangent in radians.
        module procedure :: ad_atan, ad_atan2, ad_atan2_real, real_atan2_ad
    end interface
    interface atan2
        !! Quadrant-aware inverse tangent ATAN2(y,x); not defined at the origin.
        module procedure :: ad_atan2, ad_atan2_real, real_atan2_ad
    end interface
    interface sinh
        !! Hyperbolic sine.
        module procedure :: ad_sinh
    end interface
    interface cosh
        !! Hyperbolic cosine.
        module procedure :: ad_cosh
    end interface
    interface tanh
        !! Hyperbolic tangent.
        module procedure :: ad_tanh
    end interface
    interface asinh
        !! Inverse hyperbolic sine.
        module procedure :: ad_asinh
    end interface
    interface acosh
        !! Inverse hyperbolic cosine for values >= 1; active derivative singular at 1.
        module procedure :: ad_acosh
    end interface
    interface atanh
        !! Inverse hyperbolic tangent for values strictly between -1 and 1.
        module procedure :: ad_atanh
    end interface
    interface erf
        !! Error function.
        module procedure :: ad_erf
    end interface
    interface erfc
        !! Complementary error function.
        module procedure :: ad_erfc
    end interface
    interface min
        !! Binary minimum, carrying the selected derivative; first argument wins ties.
        module procedure :: ad_min, ad_min_real, real_min_ad
    end interface
    interface max
        !! Binary maximum, carrying the selected derivative; first argument wins ties.
        module procedure :: ad_max, ad_max_real, real_max_ad
    end interface
    interface sign
        !! Magnitude of a with the sign of b; ignores the derivative of b.
        module procedure :: ad_sign, ad_sign_real, real_sign_ad
    end interface

contains
! ******************************************************************************
! ASSIGNMENT OPERATORS
! ------------------------------------------------------------------------------
elemental subroutine ad_assign_1(lhs, rhs)
    !! Copy the value and derivative of rhs to lhs.
    type(dual), intent(out) :: lhs
    !! Output dual value and derivative produced by this operation.
    type(dual), intent(in) :: rhs
    !! Input dual operand; its value and derivative participate in the documented operation.

    lhs%x = rhs%x
    lhs%dx = rhs%dx
end subroutine

elemental subroutine ad_assign_2(lhs, rhs)
    !! Assign a real64 constant rhs to lhs and clear its derivative.
    type(dual), intent(out) :: lhs
    !! Output dual value and derivative produced by this operation.
    real(real64), intent(in) :: rhs
    !! Input scalar operand or parameter used by the documented operation.
    lhs%x = rhs
    lhs%dx = 0.0d0
end subroutine

elemental subroutine ad_assign_3(lhs, rhs)
    !! Extract the value of rhs into real64 lhs, discarding the derivative.
    real(real64), intent(out) :: lhs
    !! Output scalar receiving the computed value.
    type(dual), intent(in) :: rhs
    !! Input dual operand; its value and derivative participate in the documented operation.
    lhs = rhs%x
end subroutine

elemental subroutine ad_assign_integer(lhs, rhs)
    !! Assign a default-integer constant rhs to lhs and clear its derivative.
    type(dual), intent(out) :: lhs
    !! Output dual value and derivative produced by this operation.
    integer, intent(in) :: rhs
    !! Input scalar operand or parameter used by the documented operation.
    lhs%x = real(rhs, real64)
    lhs%dx = 0.0d0
end subroutine

! ******************************************************************************
! MATHEMATICAL OPERATORS
! ------------------------------------------------------------------------------
elemental function ad_positive(a) result(c)
    !! Unary plus preserves both components.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = a
end function

! ------------------------------------------------------------------------------
elemental function ad_negative(a) result(c)
    !! Unary minus negates both components.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(-a%x, -a%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_add(a, b) result(c)
    !! Add values and directional derivatives.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x + b%x, a%dx + b%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_add_real(a, b) result(c)
    !! Add a real64 constant without changing the derivative.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x + b, a%dx)
end function

! ------------------------------------------------------------------------------
elemental function duald_ad(a, b) result(c)
    !! Add an AD operand to a real64 constant.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = b + a
end function

! ------------------------------------------------------------------------------
elemental function ad_add_integer(a, b) result(c)
    !! Add a default-integer constant using real64 arithmetic.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = a + real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental function integer_add_ad(a, b) result(c)
    !! Add an AD operand to a default-integer constant.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = real(a, real64) + b
end function

! ------------------------------------------------------------------------------
elemental function ad_subtract(a, b) result(c)
    !! Subtract values and directional derivatives.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x - b%x, a%dx - b%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_subtract_real(a, b) result(c)
    !! Subtract a real64 constant without changing the derivative.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x - b, a%dx)
end function

! ------------------------------------------------------------------------------
elemental function real_subtract_ad(a, b) result(c)
    !! Subtract an AD operand from a real64 constant, negating its derivative.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a - b%x, -b%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_subtract_integer(a, b) result(c)
    !! Subtract a default-integer constant using real64 arithmetic.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = a - real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental function integer_subtract_ad(a, b) result(c)
    !! Subtract an AD operand from a default-integer constant.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = real(a, real64) - b
end function

! ------------------------------------------------------------------------------
elemental function ad_multiply(a, b) result(c)
    !! Apply the product rule: d(a*b) = da*b + a*db.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x * b%x, a%dx * b%x + a%x * b%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_multiply_real(a, b) result(c)
    !! Scale the value and derivative by a real64 constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x * b, a%dx * b)
end function

! ------------------------------------------------------------------------------
elemental function real_multiply_ad(a, b) result(c)
    !! Multiply a real64 constant by an AD operand.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = b * a
end function

! ------------------------------------------------------------------------------
elemental function ad_multiply_integer(a, b) result(c)
    !! Scale an AD operand by a default-integer constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = a * real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental function integer_multiply_ad(a, b) result(c)
    !! Multiply a default-integer constant by an AD operand.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = real(a, real64) * b
end function

! ------------------------------------------------------------------------------
elemental function ad_divide(a, b) result(c)
    !! Apply the quotient rule as d(a/b) = (da - (a/b)*db)/b.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = a%x / b%x
    c%dx = (a%dx - c%x * b%dx) / b%x
end function

! ------------------------------------------------------------------------------
elemental function ad_divide_real(a, b) result(c)
    !! Divide the value and derivative by a nonzero real64 constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(a%x / b, a%dx / b)
end function

! ------------------------------------------------------------------------------
elemental function real_divide_ad(a, b) result(c)
    !! Divide a real64 constant by AD b, with derivative -a*db/b**2.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = a / b%x
    c%dx = -c%x * (b%dx / b%x)
end function

! ------------------------------------------------------------------------------
elemental function ad_divide_integer(a, b) result(c)
    !! Divide an AD operand by a nonzero default-integer constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = a / real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental function integer_divide_ad(a, b) result(c)
    !! Divide a default-integer constant by an AD operand.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = real(a, real64) / b
end function

! ------------------------------------------------------------------------------
elemental function ad_power(a, b) result(c)
    !! Apply d(a**b) = b*a**(b-1)*da + a**b*log(a)*db.
    !! Omit inactive terms so constant powers work at zero and negative bases.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = a ** b%x
    ! A constant exponent must not introduce log(a), which would exclude
    ! negative bases even for integral powers.
    if (b%dx /= 0.0d0) then
        if (a%x /= 0.0d0 .or. b%x <= 0.0d0) then
            c%dx = c%dx + c%x * log(a%x) * b%dx
        end if
    end if
end function

! ------------------------------------------------------------------------------
elemental function ad_power_real(a, b) result(c)
    !! Raise AD a to constant b with derivative b*a**(b-1)*da.
    !! Exponents zero and one are handled without singular intermediate powers.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = a%x ** b
    if (b == 0.0d0 .or. a%dx == 0.0d0) then
        c%dx = 0.0d0
    else if (b == 1.0d0) then
        c%dx = a%dx
    else
        c%dx = b * (a%x ** (b - 1.0d0)) * a%dx
    end if
end function

! ------------------------------------------------------------------------------
elemental function real_power_ad(a, b) result(c)
    !! Raise constant a to AD b with derivative a**b*log(a)*db.
    !! A zero base and positive exponent give zero for both components.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = a ** b%x
    if (b%dx == 0.0d0) then
        c%dx = 0.0d0
    else if (a == 0.0d0 .and. b%x > 0.0d0) then
        c%dx = 0.0d0
    else
        c%dx = c%x * log(a) * b%dx
    end if
end function

! ------------------------------------------------------------------------------
elemental function ad_power_integer(a, b) result(c)
    !! Evaluate an integer power, preserving negative-base and zero-base behavior.
    !! Negative exponents avoid b-1, which could overflow the integer range.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = a%x ** b
    if (b == 0 .or. a%dx == 0.0d0) then
        c%dx = 0.0d0
    else if (b > 0) then
        c%dx = real(b, real64) * (a%x ** (b - 1)) * a%dx
    else
        c%dx = real(b, real64) * (c%x / a%x) * a%dx
    end if
end function

! ------------------------------------------------------------------------------
elemental function integer_power_ad(a, b) result(c)
    !! Raise a default-integer constant to an AD exponent using real64 arithmetic.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = real(a, real64) ** b
end function

! ******************************************************************************
! COMPARISON OPERATORS
! ------------------------------------------------------------------------------
elemental logical function ad_equal(a, b) result(c)
    !! Compare only the values of two AD operands for equality.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a%x == b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_equal_real(a, b) result(c)
    !! Compare an AD value with a real64 constant for equality.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x == b
end function

! ------------------------------------------------------------------------------
elemental logical function real_equal_ad(a, b) result(c)
    !! Compare a real64 constant with an AD value for equality.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a == b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_equal_integer(a, b) result(c)
    !! Compare an AD value with a default integer for equality.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x == real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental logical function integer_equal_ad(a, b) result(c)
    !! Compare a default integer with an AD value for equality.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = real(a, real64) == b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_not_equal(a, b) result(c)
    !! Compare only the values of two AD operands for inequality.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a%x /= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_not_equal_real(a, b) result(c)
    !! Compare an AD value with a real64 constant for inequality.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x /= b
end function

! ------------------------------------------------------------------------------
elemental logical function real_not_equal_ad(a, b) result(c)
    !! Compare a real64 constant with an AD value for inequality.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a /= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_not_equal_integer(a, b) result(c)
    !! Compare an AD value with a default integer for inequality.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x /= real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental logical function integer_not_equal_ad(a, b) result(c)
    !! Compare a default integer with an AD value for inequality.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = real(a, real64) /= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_less(a, b) result(c)
    !! Test whether AD a's value is strictly less than AD b's.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a%x < b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_less_real(a, b) result(c)
    !! Test whether an AD value is less than a real64 constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x < b
end function

! ------------------------------------------------------------------------------
elemental logical function real_less_ad(a, b) result(c)
    !! Test whether a real64 constant is less than an AD value.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a < b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_less_integer(a, b) result(c)
    !! Test whether an AD value is less than a default integer.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x < real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental logical function integer_less_ad(a, b) result(c)
    !! Test whether a default integer is less than an AD value.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = real(a, real64) < b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_less_equal(a, b) result(c)
    !! Test whether AD a's value is less than or equal to AD b's.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a%x <= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_less_equal_real(a, b) result(c)
    !! Test whether an AD value is at most a real64 constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x <= b
end function

! ------------------------------------------------------------------------------
elemental logical function real_less_equal_ad(a, b) result(c)
    !! Test whether a real64 constant is at most an AD value.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a <= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_less_equal_integer(a, b) result(c)
    !! Test whether an AD value is at most a default integer.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x <= real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental logical function integer_less_equal_ad(a, b) result(c)
    !! Test whether a default integer is at most an AD value.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = real(a, real64) <= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_greater(a, b) result(c)
    !! Test whether AD a's value is strictly greater than AD b's.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a%x > b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_greater_real(a, b) result(c)
    !! Test whether an AD value exceeds a real64 constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x > b
end function

! ------------------------------------------------------------------------------
elemental logical function real_greater_ad(a, b) result(c)
    !! Test whether a real64 constant exceeds an AD value.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a > b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_greater_integer(a, b) result(c)
    !! Test whether an AD value exceeds a default integer.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x > real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental logical function integer_greater_ad(a, b) result(c)
    !! Test whether a default integer exceeds an AD value.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = real(a, real64) > b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_greater_equal(a, b) result(c)
    !! Test whether AD a's value is greater than or equal to AD b's.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a%x >= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_greater_equal_real(a, b) result(c)
    !! Test whether an AD value is at least a real64 constant.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x >= b
end function

! ------------------------------------------------------------------------------
elemental logical function real_greater_equal_ad(a, b) result(c)
    !! Test whether a real64 constant is at least an AD value.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = a >= b%x
end function

! ------------------------------------------------------------------------------
elemental logical function ad_greater_equal_integer(a, b) result(c)
    !! Test whether an AD value is at least a default integer.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    integer, intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    c = a%x >= real(b, real64)
end function

! ------------------------------------------------------------------------------
elemental logical function integer_greater_equal_ad(a, b) result(c)
    !! Test whether a default integer is at least an AD value.
    integer, intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    c = real(a, real64) >= b%x
end function

! ******************************************************************************
! INTRINSIC FUNCTIONS
! ------------------------------------------------------------------------------
elemental function ad_abs(a) result(c)
    !! Propagate sign(a)*da, choosing derivative zero at a = 0.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = abs(a%x)
    if (a%x == 0.0d0) then
        c%dx = 0.0d0
    else
        c%dx = sign(1.0d0, a%x) * a%dx
    end if
end function

! ------------------------------------------------------------------------------
elemental function ad_sqrt(a) result(c)
    !! Propagate da/(2*sqrt(a)); inactive seeds remain zero at a = 0.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = sqrt(a%x)
    c%dx = 0.0d0
    if (a%dx /= 0.0d0) c%dx = 0.5d0 * a%dx / c%x
end function

! ------------------------------------------------------------------------------
elemental function ad_exp(a) result(c)
    !! Propagate exp(a)*da.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = exp(a%x)
    c%dx = c%x * a%dx
end function

! ------------------------------------------------------------------------------
elemental function ad_log(a) result(c)
    !! Propagate da/a for the natural logarithm.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(log(a%x), a%dx / a%x)
end function

! ------------------------------------------------------------------------------
elemental function ad_log10(a) result(c)
    !! Propagate da/(a*log(10)) for the base-ten logarithm.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(log10(a%x), (a%dx / a%x) / log(10.0d0))
end function

! ------------------------------------------------------------------------------
elemental function ad_sin(a) result(c)
    !! Propagate cos(a)*da.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(sin(a%x), cos(a%x) * a%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_cos(a) result(c)
    !! Propagate -sin(a)*da.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(cos(a%x), -sin(a%x) * a%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_tan(a) result(c)
    !! Propagate da/cos(a)**2.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(tan(a%x), a%dx / cos(a%x)**2)
end function

! ------------------------------------------------------------------------------
elemental function ad_asin(a) result(c)
    !! Propagate da/sqrt(1-a**2), omitting the singular term for inactive seeds.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = asin(a%x)
    c%dx = 0.0d0
    if (a%dx /= 0.0d0) c%dx = a%dx / sqrt((1.0d0 - a%x) * (1.0d0 + a%x))
end function

! ------------------------------------------------------------------------------
elemental function ad_acos(a) result(c)
    !! Propagate -da/sqrt(1-a**2) using the inverse-sine derivative.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = asin(a)
    c%x = acos(a%x)
    c%dx = -c%dx
end function

! ------------------------------------------------------------------------------
elemental function ad_atan(a) result(c)
    !! Propagate da/(1+a**2), using reciprocal scaling for large a.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    real(real64) :: x_inverse
    !! Local numerical intermediate used by the derivative calculation.
    c%x = atan(a%x)
    if (abs(a%x) <= 1.0d0) then
        c%dx = a%dx / (1.0d0 + a%x * a%x)
    else
        x_inverse = safe_recip(abs(a%x))
        if (x_inverse < epsilon(1.0d0)) then
            c%dx = tiny(1.0d0)
        else
            c%dx = a%dx * x_inverse
            c%dx = c%dx * x_inverse
            c%dx = c%dx / (1.0d0 + x_inverse * x_inverse)
        end if
    end if
end function

! ------------------------------------------------------------------------------
elemental function safe_recip(x) result(r)
    !! Compute a reciprocal by scaling the fraction and exponent to reduce overflow risk.
    use ieee_arithmetic
    !! Supplies IEEE infinity construction for the zero-input case.
    real(real64), intent(in) :: x
    !! Nonzero value whose reciprocal is required.
    real(real64) :: r
    !! Computed reciprocal, or positive infinity when x is zero.
    real(real64) :: m
    !! Normalized fraction of x used in the scaled reciprocal.
    integer :: e
    !! Binary exponent of x used to restore the scale.

    if (x == 0.0d0) then
        r = ieee_value(x, ieee_positive_inf)
        return
    end if

    m = fraction(x)
    e = exponent(x)

    r = scale(1.0d0 / m, -e)
end function

! ------------------------------------------------------------------------------
elemental function ad_atan2(y, x) result(c)
    !! Propagate (x*dy-y*dx)/(x**2+y**2) using normalized coordinates.
    type(dual), intent(in) :: y, x
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    real(real64) :: scale, xscaled, yscaled
    !! Local numerical intermediate used by the derivative calculation.
    scale = max(abs(x%x), abs(y%x))
    xscaled = x%x / scale
    yscaled = y%x / scale
    c%x = atan2(y%x, x%x)
    c%dx = (xscaled * (y%dx / scale) - yscaled * (x%dx / scale)) / &
        (xscaled * xscaled + yscaled * yscaled)
end function

! ------------------------------------------------------------------------------
elemental function ad_atan2_real(y, x) result(c)
    !! Evaluate ATAN2 with AD ordinate y and constant abscissa x.
    type(dual), intent(in) :: y
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: x
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = atan2(y, dual(x, 0.0d0))
end function

! ------------------------------------------------------------------------------
elemental function real_atan2_ad(y, x) result(c)
    !! Evaluate ATAN2 with constant ordinate y and AD abscissa x.
    real(real64), intent(in) :: y
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: x
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = atan2(dual(y, 0.0d0), x)
end function

! ------------------------------------------------------------------------------
elemental function ad_sinh(a) result(c)
    !! Propagate cosh(a)*da.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(sinh(a%x), cosh(a%x) * a%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_cosh(a) result(c)
    !! Propagate sinh(a)*da.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(cosh(a%x), sinh(a%x) * a%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_tanh(a) result(c)
    !! Propagate da/cosh(a)**2 using exp(-abs(a)) to avoid hyperbolic overflow.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    real(real64) :: t, sech
    !! Local numerical intermediate used by the derivative calculation.
    t = exp(-abs(a%x))
    sech = 2.0d0 * t / (1.0d0 + t*t)
    c%x = tanh(a%x)
    c%dx = (a%dx * sech) * sech
end function

! ------------------------------------------------------------------------------
elemental function ad_asinh(a) result(c)
    !! Propagate da/sqrt(a**2+1), using HYPOT to avoid overflow.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(asinh(a%x), a%dx / hypot(a%x, 1.0d0))
end function

! ------------------------------------------------------------------------------
elemental function ad_acosh(a) result(c)
    !! Propagate da/(sqrt(a-1)*sqrt(a+1)); inactive seeds stay zero at a = 1.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c%x = acosh(a%x)
    c%dx = 0.0d0
    if (a%dx /= 0.0d0) c%dx = (a%dx / sqrt(a%x - 1.0d0)) / sqrt(a%x + 1.0d0)
end function

! ------------------------------------------------------------------------------
elemental function ad_atanh(a) result(c)
    !! Propagate da/(1-a**2) using factored divisions.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(atanh(a%x), (a%dx / (1.0d0 - a%x)) / (1.0d0 + a%x))
end function

! ------------------------------------------------------------------------------
elemental function ad_erf(a) result(c)
    !! Propagate 2/sqrt(pi)*exp(-a**2)*da.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(erf(a%x), 2.0d0 / sqrt(acos(-1.0d0)) * exp(-a%x**2) * a%dx)
end function

! ------------------------------------------------------------------------------
elemental function ad_erfc(a) result(c)
    !! Negate the error-function derivative for the complementary error function.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = erf(a)
    c%x = erfc(a%x)
    c%dx = -c%dx
end function

! ------------------------------------------------------------------------------
elemental function ad_min(a, b) result(c)
    !! Return the lower-valued operand with its derivative, choosing a on a tie.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    if (a%x <= b%x) then
        c = a
    else
        c = b
    end if
end function

! ------------------------------------------------------------------------------
elemental function ad_min_real(a, b) result(c)
    !! Select the minimum of AD a and constant b; selecting b clears the derivative.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = min(a, dual(b, 0.0d0))
end function

! ------------------------------------------------------------------------------
elemental function real_min_ad(a, b) result(c)
    !! Select the minimum of constant a and AD b; selecting a clears the derivative.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = min(dual(a, 0.0d0), b)
end function

! ------------------------------------------------------------------------------
elemental function ad_max(a, b) result(c)
    !! Return the higher-valued operand with its derivative, choosing a on a tie.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    if (a%x >= b%x) then
        c = a
    else
        c = b
    end if
end function

! ------------------------------------------------------------------------------
elemental function ad_max_real(a, b) result(c)
    !! Select the maximum of AD a and constant b; selecting b clears the derivative.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = max(a, dual(b, 0.0d0))
end function

! ------------------------------------------------------------------------------
elemental function real_max_ad(a, b) result(c)
    !! Select the maximum of constant a and AD b; selecting a clears the derivative.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = max(dual(a, 0.0d0), b)
end function

! ------------------------------------------------------------------------------
elemental function ad_sign(a, b) result(c)
    !! Propagate sign(1,b)*d(abs(a)), treating the sign selection as constant.
    type(dual), intent(in) :: a, b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = abs(a)
    c%x = sign(a%x, b%x)
    c%dx = sign(1.0d0, b%x) * c%dx
end function

! ------------------------------------------------------------------------------
elemental function ad_sign_real(a, b) result(c)
    !! Apply a constant sign to AD a and its absolute-value derivative.
    type(dual), intent(in) :: a
    !! Input dual operand; its value and derivative participate in the documented operation.
    real(real64), intent(in) :: b
    !! Input scalar operand or parameter used by the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = sign(a, dual(b, 0.0d0))
end function

! ------------------------------------------------------------------------------
elemental function real_sign_ad(a, b) result(c)
    !! Apply AD b's sign to a constant magnitude, returning derivative zero.
    real(real64), intent(in) :: a
    !! Input scalar operand or parameter used by the documented operation.
    type(dual), intent(in) :: b
    !! Input dual operand; its value and derivative participate in the documented operation.
    type(dual) :: c
    !! Computed dual-number result containing a value and propagated derivative.
    c = dual(sign(a, b%x), 0.0d0)
end function

! ------------------------------------------------------------------------------
end module
