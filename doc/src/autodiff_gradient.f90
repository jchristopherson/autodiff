module autodiff_gradient
    !! Forward-mode evaluation of scalar-function gradients.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    use autodiff_workspace
    use autodiff_hyperdual_numbers, only: hyperdual, assignment(=)
    use autodiff_hyperdual_workspace, only: hyperdual_workspace
    implicit none

    private :: gradient_dual, gradient_hyperdual
    private :: hyperdual, hyperdual_workspace

    interface gradient
        !! Compute a scalar function's gradient with a dual or hyper-dual callback.
        module procedure :: gradient_dual, gradient_hyperdual
    end interface

contains

    subroutine gradient_dual(fcn, x, g, work, args)
        !! Compute a scalar function's gradient by seeding one input direction 
        !! at a time.
        procedure(scalar_function) :: fcn
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

    subroutine gradient_hyperdual(fcn, x, g, work, args)
        !! Compute a scalar function's gradient using the first hyper-dual
        !! derivative direction, leaving the other seeds zero.
        procedure(hyperdual_function) :: fcn
            !! User function whose result carries the value and derivatives.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which the gradient is evaluated.
        real(real64), intent(out), dimension(:) :: g
            !! Gradient output, with one component per input.
        class(hyperdual_workspace), intent(inout), target, optional :: work
            !! Reusable input buffer, grown if needed.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the function.

        type(hyperdual), allocatable, target, dimension(:) :: x_buffer
        type(hyperdual), pointer, dimension(:) :: xd
        type(hyperdual) :: fd
        integer :: i, n

        n = size(x)
        if (present(work)) then
            if (.not.work%validate_size(n)) call work%resize(n)
            xd => work%xd_pool(1:n, 1)
        else
            allocate(x_buffer(n))
            xd => x_buffer
        end if

        do i = 1, n
            xd(i) = x(i)
        end do

        do i = 1, n
            xd(i)%dx = 1.0_real64
            fd = fcn(xd, args = args)
            g(i) = fd%dx
            xd(i)%dx = 0.0_real64
        end do
    end subroutine

end module