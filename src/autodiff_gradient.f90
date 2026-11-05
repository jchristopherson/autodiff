module autodiff_gradient
    !! Forward-mode evaluation of scalar-function gradients.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    implicit none

contains

    subroutine gradient(fcn, x, g, args)
        !! Compute a scalar function's gradient by seeding one input direction 
        !! at a time.
        procedure(scalar_function), pointer, intent(in) :: fcn
            !! User function whose result carries the value and seeded derivative.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which the gradient is evaluated.
        real(real64), intent(out), dimension(:) :: g
            !! Gradient output, with one component per input; 
            !! \(g_i = \partial f/\partial x_i\).
        class(*), intent(inout), optional :: args
            !! Optional caller-owned context passed to the procedure, which 
            !! may update it.

        type(dual), allocatable, dimension(:) :: xd
        type(dual) :: fd
        integer :: i, j, n

        n = size(x)
        allocate(xd(n))

        do i = 1, n
            do j = 1, n
                xd(j)%x  = x(j)
                xd(j)%dx = 0.0_real64
            end do

            xd(i)%dx = 1.0_real64
            fd = fcn(xd, args = args)
            g(i) = fd%dx
        end do
    end subroutine

end module