module differentiation_test_support
    use autodiff_dual_numbers, only: dual, cos, operator(+), operator(*)
    implicit none

    type :: call_counter
        integer, pointer :: calls => null()
    end type

contains

    subroutine vector_test_function(xd, y, args)
        type(dual), intent(in) :: xd(:)
        type(dual), intent(out) :: y(:)
        class(*), intent(in), optional :: args

        if (present(args)) then
            select type (args)
            type is (call_counter)
                !$omp atomic
                args%calls = args%calls + 1
            class default
                error stop "Unexpected vector-function context type."
            end select
        end if
        y(1) = xd(1) * xd(2)
        y(2) = xd(1) + cos(xd(2))
    end subroutine

end module

program test_differentiation
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_dual_numbers, only: dual, operator(+), operator(*), operator(**), sin, cos
    use differentiation_test_support, only: call_counter, vector_test_function
    use autodiff_gradient, only: gradient
    use autodiff_jacobian, only: jacobian
    use autodiff_jvp, only: jvp
    use autodiff_interfaces, only: scalar_function, vector_function
    use autodiff_workspace, only: workspace
    implicit none

    integer, target :: total_calls
    type(call_counter) :: counter
    procedure(scalar_function), pointer :: scalar_proc
    procedure(vector_function), pointer :: vector_proc
    real(real64), parameter :: x(2) = [2.0_real64, 0.5_real64]
    real(real64), parameter :: v(2) = [3.0_real64, -2.0_real64]
    real(real64) :: g(2), j(2,2), jv_result(2)
    real(real64) :: j_ws(2,2), jv_ws(2), j_par(2,2)
    type(workspace) :: work, work_new, work_large, work_small, work_jvp
    type(workspace) :: work_grad, work_par
    integer :: checks

    checks = 0
    total_calls = 0
    counter%calls => total_calls
    scalar_proc => scalar_test_function
    vector_proc => vector_test_function

    call gradient(scalar_proc, x, g, args = counter)
    call check_close("gradient x1", g(1), 2.0_real64 * x(1) * x(2))
    call check_close("gradient x2", g(2), x(1)**2 + cos(x(2)))
    call check_true("gradient forwards optional context", counter%calls == 2)

    call gradient(scalar_proc, x, g)
    call check_true("gradient supports absent optional context", counter%calls == 2)

    call jacobian(vector_proc, x, j, args = counter)
    call check_close("Jacobian (1,1)", j(1,1), x(2))
    call check_close("Jacobian (1,2)", j(1,2), x(1))
    call check_close("Jacobian (2,1)", j(2,1), 1.0_real64)
    call check_close("Jacobian (2,2)", j(2,2), -sin(x(2)))
    call check_true("Jacobian forwards optional context", counter%calls == 4)

    call jacobian(vector_proc, x, j)
    call check_true("Jacobian supports absent optional context", counter%calls == 4)

    call jvp(vector_proc, x, v, jv_result, args = counter)
    call check_close("JVP output 1", jv_result(1), x(2) * v(1) + x(1) * v(2))
    call check_close("JVP output 2", jv_result(2), v(1) - sin(x(2)) * v(2))
    call check_true("JVP forwards optional context", counter%calls == 5)

    call jvp(vector_proc, x, v, jv_result)
    call check_true("JVP supports absent optional context", counter%calls == 5)

    call work%resize(2, 3)
    call check_true("workspace allocates requested buffers", &
        allocated(work%xd) .and. allocated(work%fd))
    call check_true("workspace input size", size(work%xd) == 2)
    call check_true("workspace output size", size(work%fd) == 3)
    call work%resize(2, 3)
    call check_true("workspace preserves matching buffers", &
        size(work%xd) == 2 .and. size(work%fd) == 3)
    call work%resize(4, 3)
    call check_true("workspace resizes inputs independently", size(work%xd) == 4 .and. size(work%fd) == 3)
    call work%resize(4, 1)
    call check_true("workspace never shrinks outputs", size(work%xd) == 4 .and. size(work%fd) == 3)

    call work_large%resize(5, 4)
    call check_true("workspace validates larger buffers", work_large%validate_size(2, 2))
    call check_true("workspace rejects too-small buffers", .not. work_large%validate_size(6, 2))

    call jacobian(vector_proc, x, j_ws, work_new, args = counter)
    call check_true("Jacobian with workspace forwards context", counter%calls == 7)
    call check_close("Jacobian with workspace (1,1)", j_ws(1,1), x(2))
    call check_close("Jacobian with workspace (1,2)", j_ws(1,2), x(1))
    call check_close("Jacobian with workspace (2,1)", j_ws(2,1), 1.0_real64)
    call check_close("Jacobian with workspace (2,2)", j_ws(2,2), -sin(x(2)))
    call check_true("Jacobian allocates empty workspace", &
        size(work_new%xd) == 2 .and. size(work_new%fd) == 2 .and. &
        size(work_new%xd_pool, 2) == 1 .and. size(work_new%fd_pool, 2) == 1)
    call check_close("Jacobian workspace input seed 1 cleared", work_new%xd_pool(1,1)%dx, 0.0_real64)
    call check_close("Jacobian workspace input seed 2 cleared", work_new%xd_pool(2,1)%dx, 0.0_real64)
    call check_close("Jacobian workspace output value", work_new%fd_pool(1,1)%x, x(1) * x(2))
    call check_close("Jacobian workspace output derivative", work_new%fd_pool(2,1)%dx, -sin(x(2)))

    call jacobian(vector_proc, x, j_ws, work_new)
    call check_true("Jacobian without context forwards no calls", counter%calls == 7)
    call check_close("Jacobian reused workspace (2,1)", j_ws(2,1), 1.0_real64)

    call jacobian(vector_proc, x, j_ws, work_large)
    call check_true("Jacobian preserves oversized workspace", &
        size(work_large%xd) == 5 .and. size(work_large%fd) == 4)
    call check_close("Jacobian oversized workspace (1,2)", j_ws(1,2), x(1))
    call check_close("Jacobian oversized workspace (2,2)", j_ws(2,2), -sin(x(2)))

    call work_small%resize(1, 1)
    call jacobian(vector_proc, x, j_ws, work_small)
    call check_true("Jacobian grows undersized workspace", &
        size(work_small%xd) == 2 .and. size(work_small%fd) == 2)
    call check_close("Jacobian undersized workspace (1,1)", j_ws(1,1), x(2))

    call jvp(vector_proc, x, v, jv_ws, work_jvp, args = counter)
    call check_true("JVP with workspace forwards context", counter%calls == 8)
    call check_close("JVP with workspace output 1", jv_ws(1), x(2) * v(1) + x(1) * v(2))
    call check_close("JVP with workspace output 2", jv_ws(2), v(1) - sin(x(2)) * v(2))
    call check_true("JVP allocates empty workspace", &
        size(work_jvp%xd) == 2 .and. size(work_jvp%fd) == 2)
    call check_close("JVP workspace input seed 1", work_jvp%xd(1)%dx, v(1))
    call check_close("JVP workspace input seed 2", work_jvp%xd(2)%dx, v(2))
    call check_close("JVP workspace output derivative", work_jvp%fd(2)%dx, jv_ws(2))

    call jvp(vector_proc, x, v, jv_ws, work_large)
    call check_true("JVP preserves oversized workspace", &
        size(work_large%xd) == 5 .and. size(work_large%fd) == 4)
    call check_close("JVP oversized workspace output 1", jv_ws(1), x(2) * v(1) + x(1) * v(2))
    call check_close("JVP oversized workspace output 2", jv_ws(2), v(1) - sin(x(2)) * v(2))
    call check_close("JVP oversized workspace input seed 2", work_large%xd(2)%dx, v(2))

    call work_small%resize(1, 1)
    call jvp(vector_proc, x, v, jv_ws, work_small)
    call check_true("JVP grows undersized workspace", &
        size(work_small%xd) == 2 .and. size(work_small%fd) == 2)
    call check_close("JVP undersized workspace output 1", jv_ws(1), x(2) * v(1) + x(1) * v(2))
    call check_close("JVP undersized workspace output 2", jv_ws(2), v(1) - sin(x(2)) * v(2))

    call gradient(scalar_proc, x, g, work = work_grad, args = counter)
    call check_close("gradient with workspace x1", g(1), 2.0_real64 * x(1) * x(2))
    call check_close("gradient with workspace x2", g(2), x(1)**2 + cos(x(2)))
    call check_true("gradient with workspace forwards context", counter%calls == 10)
    call check_true("gradient allocates input workspace", size(work_grad%xd) == 2)
    call check_close("gradient workspace seed cleared", work_grad%xd(2)%dx, 0.0_real64)

    ! call work_par%resize(2, 2, 2)
    ! call check_true("parallel workspace allocates thread columns", &
    !     size(work_par%xd_pool, 2) == 2 .and. size(work_par%fd_pool, 2) == 2)
    ! call jacobian(vector_proc, x, j_par, work_par, nthreads = 2)
    ! call check_close("parallel Jacobian (1,1)", j_par(1,1), x(2))
    ! call check_close("parallel Jacobian (1,2)", j_par(1,2), x(1))
    ! call check_close("parallel Jacobian (2,1)", j_par(2,1), 1.0_real64)
    ! call check_close("parallel Jacobian (2,2)", j_par(2,2), -sin(x(2)))

    ! call jacobian(vector_proc, x, j_par, nthreads = 2)
    ! call check_close("parallel Jacobian without workspace (1,2)", j_par(1,2), x(1))
    ! call check_close("parallel Jacobian without workspace (2,2)", j_par(2,2), -sin(x(2)))

    write(output_unit, '(a,i0,a)') "Passed ", checks, " differentiation checks."

contains

    function scalar_test_function(xd, args) result(y)
        type(dual), intent(in) :: xd(:)
        class(*), intent(in), optional :: args
        type(dual) :: y

        if (present(args)) then
            select type (args)
            type is (call_counter)
                !$omp atomic
                args%calls = args%calls + 1
            class default
                error stop "Unexpected scalar-function context type."
            end select
        end if
        y = xd(1)**2 * xd(2) + sin(xd(2))
    end function

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
