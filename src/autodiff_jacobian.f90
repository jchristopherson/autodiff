module autodiff_jacobian
    !! Forward-mode evaluation of vector-function Jacobian matrices.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    use autodiff_workspace
    implicit none

contains

    subroutine jacobian(fcn, x, J, work, args)
        !! Compute a vector function's Jacobian by seeding each input basis 
        !! direction.
        procedure(vector_function), pointer, intent(in) :: fcn
            !! User procedure that evaluates the vector-valued function.
        real(real64), intent(in), dimension(:) :: x
            !! Input point at which the Jacobian is evaluated.
        real(real64), intent(out), dimension(:,:) :: J
            !! Jacobian output with shape (number of outputs, number of inputs);
            !! \(J_{ik} = \partial f_i/\partial x_k\).
        class(workspace), intent(inout), target, optional :: work
            !! An optional workspace to use in order to avoid internal memory
            !! allocations.  If not sized properly, this routine will reallocate
            !! the workspace to suit its needs.
        class(*), intent(inout), optional :: args
            !! Optional caller-owned context passed to the procedure, which 
            !! may update it.

        type(dual), allocatable, target, dimension(:) :: x_buffer, f_buffer
        type(dual), pointer, dimension(:) :: xd, fd
        integer :: i, k, m, n

        n = size(x)
        m = size(J,1)

        if (present(work)) then
            if (.not.work%validate_size(n, m)) call work%resize(n, m)
            xd => work%xd(1:n)
            fd => work%fd(1:m)
        else
            allocate(x_buffer(n), f_buffer(m))
            xd => x_buffer
            fd => f_buffer
        end if

        do k = 1, n
            do i = 1, n
                xd(i)%x = x(i)
                xd(i)%dx = 0.0_real64
            end do

            xd(k)%dx = 1.0_real64

            call fcn(xd, fd, args = args)

            do i = 1,m
                J(i,k) = fd(i)%dx
            end do
        end do

    end subroutine

end module