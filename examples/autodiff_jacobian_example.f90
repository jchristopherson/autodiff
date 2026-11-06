program example
    use iso_fortran_env
    use autodiff
    implicit none

    real(real64) :: x(2), J(2,2), J_analytical(2,2), Delta(2,2)

    ! Define the x and y values at which the Jacobian will be evaluated
    call random_number(x)

    ! Define the analytical Jacobian
    J_analytical = reshape([2.0d0 * x(1), 1.0d0, 2.0d0 * x(2), -1.0d0], [2, 2])

    ! The automatic differentiation approach
    call jacobian(vec_fcn, x, J)

    ! Compute any differences
    Delta = J - J_analytical

    ! Display the results
    print "(A)", "ANALYTICAL:"
    call write_matrix(J_analytical)

    print "(A)", new_line('a') // "AUTO DIFF:"
    call write_matrix(J)

    print "(A)", new_line('a') // "DIFFERENCES:"
    call write_matrix(Delta)
    

contains
    subroutine vec_fcn(u, f, args)
        ! Defines the functions:
        ! F1 = x**2 + y**2 - 1
        ! F2 = x - y
        type(dual), intent(in), dimension(:) :: u
        type(dual), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args

        f(1) = u(1)**2 + u(2)**2 - 1.0d0
        f(2) = u(1) - u(2)
    end subroutine

    subroutine write_matrix(z)
        real(real64), intent(in), dimension(:,:) :: z
        print "(E0.3, A, E0.3)", z(1,1), achar(9), z(1,2)
        print "(E0.3, A, E0.3)", z(2,1), achar(9), z(2,2)
    end subroutine
end program