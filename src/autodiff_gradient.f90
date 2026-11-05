module autodiff_gradient
    !! Forward-mode evaluation of scalar-function gradients.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    use autodiff_workspace
    implicit none

contains

    subroutine gradient(fcn, x, g, work, args)
        !! Compute a scalar function's gradient by seeding one input direction 
        !! at a time.
        procedure(scalar_function), pointer, intent(in) :: fcn
            !! User function whose result carries the value and seeded derivative.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which the gradient is evaluated.
        real(real64), intent(out), dimension(:) :: g
            !! Gradient output, with one component per input; 
            !! \(g_i = \partial f/\partial x_i\).
        class(workspace), intent(inout), target, optional :: work
            !! An optional workspace to use in order to avoid internal memory
            !! allocations.  If not sized properly, this routine will grow
            !! the workspace to suit its needs.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the procedure.

        type(dual), allocatable, target, dimension(:) :: x_buffer
        type(dual), pointer, dimension(:) :: xd
        type(dual) :: fd
        integer :: i, j, n

        n = size(x)

        if (present(work)) then
            if (.not.work%validate_size(n, 1)) call work%resize(n, 1)
            xd => work%xd(1:n)
        else
            allocate(x_buffer(n))
            xd => x_buffer
        end if

        do j = 1, n
            xd(j)%x  = x(j)
            xd(j)%dx = 0.0d0
        end do

        do i = 1, n
            xd(i)%dx = 1.0d0
            fd = fcn(xd, args = args)
            g(i) = fd%dx
            xd(i)%dx = 0.0d0
        end do
    end subroutine

end module