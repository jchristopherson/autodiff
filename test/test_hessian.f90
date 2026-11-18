program test_hessian
    use iso_fortran_env, only: real64, error_unit, output_unit
    use autodiff_hyperdual_numbers
    use autodiff_interfaces, only: hyperdual_function
    use autodiff_hessian, only: hessian, hvp
    use autodiff_hyperdual_workspace, only: hyperdual_workspace
    implicit none

    type :: call_counter
        integer, pointer :: calls => null()
    end type

    real(real64), parameter :: tol = 1.0d-10
    integer :: failures
    integer, target :: total_calls
    real(real64) :: x(2), H(2,2), g(2), f, v(2), Hv(2)
    real(real64) :: x3(3), H3(3,3), g3(3), f3, Hp(3,3), gp(3), fp
    real(real64) :: x0(0), H0(0,0)
    type(call_counter) :: counter
    type(hyperdual_workspace) :: work, work_par, work_large, work_small

    failures = 0
    total_calls = 0
    counter%calls => total_calls

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

    ! Hessian-vector product, checked against the analytic mixed function Hessian
    x = [1.3_real64, 0.7_real64]
    v = [0.4_real64, -0.9_real64]
    call hvp(mixed_function, x, v, Hv)
    H(1,1) = 2.0_real64 * x(2) + x(2)**2 * exp(x(1) * x(2))
    H(1,2) = 2.0_real64 * x(1) + exp(x(1) * x(2)) + x(1) * x(2) * exp(x(1) * x(2))
    H(2,2) = -sin(x(2)) + x(1)**2 * exp(x(1) * x(2))
    call check_close("hvp component 1", Hv(1), H(1,1) * v(1) + H(1,2) * v(2))
    call check_close("hvp component 2", Hv(2), H(1,2) * v(1) + H(2,2) * v(2))

    ! Value and gradient come from the diagonal evaluations: n(n+1)/2 calls
    x3 = [0.8_real64, -0.4_real64, 1.7_real64]
    call hessian(counted_function, x3, H3, g3, f3, args = counter)
    call check_true("hessian evaluation count", counter%calls == 6)
    call check_close("three-variable value", f3, x3(1)**2 * x3(2) * x3(3) + sin(x3(1) * x3(3)))
    call check_close("three-variable g3", g3(3), x3(1)**2 * x3(2) + x3(1) * cos(x3(1) * x3(3)))
    call check_close("three-variable H13", H3(1,3), 2.0_real64 * x3(1) * x3(2) &
        + cos(x3(1) * x3(3)) - x3(1) * x3(3) * sin(x3(1) * x3(3)))
    call check_close("three-variable H31", H3(3,1), H3(1,3))

    ! Workspace allocation, seed cleanup, reuse, and growth
    call hessian(counted_function, x3, Hp, gp, fp, work = work)
    call check_true("hessian allocates empty workspace", &
        size(work%xd_pool, 1) == 3 .and. size(work%xd_pool, 2) == 1)
    call check_true("hessian workspace seeds cleared", &
        all(work%xd_pool(:,1)%dx == 0.0_real64) .and. all(work%xd_pool(:,1)%dy == 0.0_real64))
    call check_matrix("hessian with workspace", Hp, H3)
    call check_close("hessian with workspace value", fp, f3)
    call check_close("hessian with workspace g2", gp(2), g3(2))
    call hessian(counted_function, x3, Hp, work = work)
    call check_matrix("hessian with reused workspace", Hp, H3)

    call work_large%resize(5, 3)
    call hessian(counted_function, x3, Hp, work = work_large)
    call check_true("hessian preserves oversized workspace", &
        size(work_large%xd_pool, 1) == 5 .and. size(work_large%xd_pool, 2) == 3)
    call check_matrix("hessian with oversized workspace", Hp, H3)

    call work_small%resize(1)
    call hessian(counted_function, x3, Hp, work = work_small, nthreads = 2)
    call check_true("hessian grows undersized workspace", &
        size(work_small%xd_pool, 1) >= 3 .and. size(work_small%xd_pool, 2) >= 2)
    call check_matrix("hessian with grown workspace", Hp, H3)

    ! Threaded evaluation must match the serial result
    counter%calls = 0
    call work_par%resize(3, 2)
    call hessian(counted_function, x3, Hp, gp, fp, work_par, args = counter, nthreads = 2)
    call check_true("parallel hessian evaluation count", counter%calls == 6)
    call check_matrix("parallel hessian", Hp, H3)
    call check_close("parallel hessian value", fp, f3)
    call check_close("parallel hessian g1", gp(1), g3(1))
    call check_close("parallel hessian g3", gp(3), g3(3))

    call hessian(counted_function, x3, Hp, gp, nthreads = 4)
    call check_matrix("parallel hessian without workspace", Hp, H3)
    call check_close("parallel hessian without workspace g2", gp(2), g3(2))

    ! Hessian-vector product with a workspace
    call hvp(mixed_function, x, v, Hv, work = work)
    call check_close("hvp with workspace component 1", Hv(1), H(1,1) * v(1) + H(1,2) * v(2))
    call check_close("hvp with workspace component 2", Hv(2), H(1,2) * v(1) + H(2,2) * v(2))

    ! An empty input still yields the function value
    counter%calls = 0
    call hessian(counted_function, x0, H0, f = f, args = counter)
    call check_true("empty-input hessian evaluates once", counter%calls == 1)
    call check_close("empty-input hessian value", f, 0.0_real64)

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

    function counted_function(u, args) result(fv)
        type(hyperdual), intent(in), dimension(:) :: u
        class(*), intent(in), optional :: args
        type(hyperdual) :: fv

        if (present(args)) then
            select type (args)
            type is (call_counter)
                !$omp atomic
                args%calls = args%calls + 1
            class default
                error stop "Unexpected hessian context type."
            end select
        end if
        if (size(u) == 0) then
            fv = 0.0_real64
        else
            fv = u(1)**2 * u(2) * u(3) + sin(u(1) * u(3))
        end if
    end function

    subroutine check_true(label, condition)
        character(len=*), intent(in) :: label
        logical, intent(in) :: condition
        if (.not. condition) then
            failures = failures + 1
            write(error_unit, '(A,A)') "FAILED: ", label
        end if
    end subroutine

    subroutine check_matrix(label, actual, expected)
        character(len=*), intent(in) :: label
        real(real64), intent(in) :: actual(:,:), expected(:,:)
        integer :: i, j
        do j = 1, size(expected, 2)
            do i = 1, size(expected, 1)
                call check_close(label, actual(i,j), expected(i,j))
            end do
        end do
    end subroutine

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
