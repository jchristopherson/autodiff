program test_differentiation
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_dual_numbers, only: dual, operator(+), operator(*), operator(**), sin, cos
    use autodiff_gradient, only: gradient
    use autodiff_jacobian, only: jacobian
    use autodiff_jvp, only: jvp
    use autodiff_interfaces, only: scalar_function, vector_function
    use autodiff_workspace, only: ad_work, resize
    implicit none

    type :: call_counter
        integer :: calls = 0
    end type

    type(call_counter) :: counter
    procedure(scalar_function), pointer :: scalar_proc
    procedure(vector_function), pointer :: vector_proc
    real(real64), parameter :: x(2) = [2.0_real64, 0.5_real64]
    real(real64), parameter :: v(2) = [3.0_real64, -2.0_real64]
    real(real64) :: g(2), j(2,2), jv_result(2)
    type(ad_work) :: work
    integer :: checks

    checks = 0
    scalar_proc => scalar_test_function
    vector_proc => vector_test_function

    call gradient(scalar_proc, x, g, counter)
    call check_close("gradient x1", g(1), 2.0_real64 * x(1) * x(2))
    call check_close("gradient x2", g(2), x(1)**2 + cos(x(2)))
    call check_true("gradient forwards optional context", counter%calls == 2)

    call gradient(scalar_proc, x, g)
    call check_true("gradient supports absent optional context", counter%calls == 2)

    call jacobian(vector_proc, x, j, counter)
    call check_close("Jacobian (1,1)", j(1,1), x(2))
    call check_close("Jacobian (1,2)", j(1,2), x(1))
    call check_close("Jacobian (2,1)", j(2,1), 1.0_real64)
    call check_close("Jacobian (2,2)", j(2,2), -sin(x(2)))
    call check_true("Jacobian forwards optional context", counter%calls == 4)

    call jacobian(vector_proc, x, j)
    call check_true("Jacobian supports absent optional context", counter%calls == 4)

    call jvp(vector_proc, x, v, jv_result, counter)
    call check_close("JVP output 1", jv_result(1), x(2) * v(1) + x(1) * v(2))
    call check_close("JVP output 2", jv_result(2), v(1) - sin(x(2)) * v(2))
    call check_true("JVP forwards optional context", counter%calls == 5)

    call jvp(vector_proc, x, v, jv_result)
    call check_true("JVP supports absent optional context", counter%calls == 5)

    call resize(work, 2, 3)
    call check_true("workspace allocates requested buffers", &
        allocated(work%xd) .and. allocated(work%fd))
    call check_true("workspace input size", size(work%xd) == 2)
    call check_true("workspace output size", size(work%fd) == 3)
    call resize(work, 2, 3)
    call check_true("workspace preserves matching buffers", &
        size(work%xd) == 2 .and. size(work%fd) == 3)
    call resize(work, 4, 3)
    call check_true("workspace resizes inputs independently", size(work%xd) == 4 .and. size(work%fd) == 3)
    call resize(work, 4, 1)
    call check_true("workspace resizes outputs independently", size(work%xd) == 4 .and. size(work%fd) == 1)

    write(output_unit, '(a,i0,a)') "Passed ", checks, " differentiation checks."

contains

    function scalar_test_function(xd, args) result(y)
        type(dual), intent(in) :: xd(:)
        class(*), intent(inout), optional :: args
        type(dual) :: y

        if (present(args)) then
            select type (args)
            type is (call_counter)
                args%calls = args%calls + 1
            class default
                error stop "Unexpected scalar-function context type."
            end select
        end if
        y = xd(1)**2 * xd(2) + sin(xd(2))
    end function

    subroutine vector_test_function(xd, y, args)
        type(dual), intent(in) :: xd(:)
        type(dual), intent(out) :: y(:)
        class(*), intent(inout), optional :: args

        if (present(args)) then
            select type (args)
            type is (call_counter)
                args%calls = args%calls + 1
            class default
                error stop "Unexpected vector-function context type."
            end select
        end if
        y(1) = xd(1) * xd(2)
        y(2) = xd(1) + cos(xd(2))
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
