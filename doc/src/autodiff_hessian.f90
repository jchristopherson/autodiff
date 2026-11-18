module autodiff_hessian
    !! Second-order forward-mode evaluation of scalar-function gradients and Hessians.
    use iso_fortran_env
    use autodiff_hyperdual_numbers
    use autodiff_interfaces
    implicit none

contains

    subroutine hessian(fcn, x, H, g, f, args)
        !! Compute a scalar function's Hessian, and optionally its gradient and
        !! value, by seeding hyper-dual inputs.  The Hessian is symmetric and
        !! requires \(n(n+1)/2\) evaluations of `fcn`, plus \(n\) more when the
        !! gradient is requested.
        procedure(hyperdual_function) :: fcn
            !! User function whose result carries the value and derivatives.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which the Hessian is evaluated.
        real(real64), intent(out), dimension(:,:) :: H
            !! Hessian output, \(H_{ij} = \partial^2 f / \partial x_i \partial x_j\).
        real(real64), intent(out), optional, dimension(:) :: g
            !! Optional gradient output, \(g_i = \partial f/\partial x_i\).
        real(real64), intent(out), optional :: f
            !! Optional function value at `x`.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the function.

        type(hyperdual), allocatable, dimension(:) :: xd
        type(hyperdual) :: fd
        integer :: i, j, n

        n = size(x)
        allocate(xd(n))
        do j = 1, n
            xd(j) = x(j)
        end do

        do i = 1, n
            xd(i)%dx = 1.0_real64
            fd = fcn(xd, args = args)
            if (i == 1 .and. present(f)) f = fd%x
            if (present(g)) g(i) = fd%dx
            xd(i)%dx = 0.0_real64

            do j = i, n
                xd(i)%dx = 1.0_real64
                xd(j)%dy = 1.0_real64
                fd = fcn(xd, args = args)
                H(i, j) = fd%dxy
                H(j, i) = fd%dxy
                xd(i)%dx = 0.0_real64
                xd(j)%dy = 0.0_real64
            end do
        end do
    end subroutine

end module
