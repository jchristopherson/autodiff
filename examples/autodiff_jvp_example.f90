program example
    use iso_fortran_env
    use autodiff
    implicit none

    ! Model Parameters
    real(real64), parameter :: a = 1.5d0
    real(real64), parameter :: b = 2.3d0
    real(real64), parameter :: c = 0.1d0
    real(real64), parameter :: d = 5.7d0

    ! Variables
    integer :: i
    real(real64) :: x(2), v(2), Jv(2), J(2,2), Ans(2), delta(2)

    ! Define x and v
    call random_number(x)
    call random_number(v)

    ! Calculate J * v directly, without forming the Jacobian J
    call jvp(vec_fcn, x, v, Jv)

    ! For comparison, form J and multiply by v directly
    call jacobian(vec_fcn, x, J)
    Ans = matmul(J, v)

    ! Compute the differences
    delta = Jv - Ans

    ! Display the results
    print "(A)", "JVP:"
    print "(E0.3)", (Jv(i), i = 1, size(Jv))

    print "(A)", new_line('a') // "FULLY FORMED J:"
    print "(E0.3)", (Ans(i), i = 1, size(Ans))

    print "(A)", new_line('a') // "DIFFERENCES:"
    print "(E0.3)", (delta(i), i = 1, size(delta))

contains
    subroutine vec_fcn(u, f, args)
        ! Consider the ODE's:
        ! x' = a * x - b x * y
        ! y' = c * x * y - d * y
        type(dual), intent(in), dimension(:) :: u
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args

        f(1) = a * u(1) - b * u(1) * u(2)
        f(2) = c * u(1) * u(2) - d * u(2)
    end subroutine
end program