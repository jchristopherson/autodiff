program example
    use iso_fortran_env
    use autodiff
    implicit none

    integer :: i
    real(real64) :: x(2), g(2), g_analytical(2), delta(2)

    ! Define the x and y values at which to evaluate the gradient vector
    call random_number(x)

    ! Define the analytical gradient
    g_analytical = [ &
        -4.0d2 * x(1) * (x(2) - x(1)**2) - 2.0d0 * (1.0d0 - x(1)), &
        2.0d2 * (x(2) - x(1)**2) &
    ]

    ! The automatic differentiation approach
    call gradient(sclr_fcn, x, g)

    ! Compute any differences
    delta = g - g_analytical

    ! Display the results
    print "(A)", "ANALYTICAL:"
    print "(E0.3)", (g_analytical(i), i = 1, size(g_analytical))

    print "(A)", new_line('a') // "AUTO DIFF:"
    print "(E0.3)", (g(i), i = 1, size(g))

    print "(A)", new_line('a') // "DIFFERENCES:"
    print "(E0.3)", (delta(i), i = 1, size(delta))
contains
    function sclr_fcn(u, args) result(f)
        ! Rosenbrock's Function:
        ! f(x,y) = (1 - x)**2 + 100 * (y - x**2)**2
        type(dual), intent(in), dimension(:) :: u
        class(*), intent(in), optional :: args
        type(dual) :: f
        f = (1.0d0 - u(1))**2 + 1.0d2 * (u(2) - u(1)**2)**2
    end function
end program