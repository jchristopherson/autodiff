module autodiff_jvp
    !! Forward-mode evaluation of Jacobian-vector products.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    use autodiff_workspace

contains

    subroutine jvp(fcn, x, v, Jv, work, args)
        !! Evaluate the directional derivative of a vector function without 
        !! forming its Jacobian.
        procedure(vector_function), pointer, intent(in) :: fcn
            !! User procedure that evaluates the vector-valued function.
        real(real64), intent(in) :: x(:)
            !! Input point at which the directional derivative is evaluated.
        real(real64), intent(in), dimension(:) :: v
            !! Direction vector, with one component per input; its entries seed 
            !! each dual input.
        real(real64), intent(out), dimension(:) :: Jv
            !! Product of the Jacobian and v, satisfying \(Jv = J_f(x)v\).
        class(workspace), intent(inout), target, optional :: work
            !! An optional workspace to use in order to avoid internal memory
            !! allocations.  If not sized properly, this routine will reallocate
            !! the workspace to suit its needs.
        class(*), intent(inout), optional :: args
            !! Optional caller-owned context passed to the procedure, which 
            !! may update it.

        type(dual), allocatable, target, dimension(:) :: x_buffer, f_buffer
        type(dual), pointer, dimension(:) :: xd, fd
        integer :: i, n, m

        n = size(x)
        m = size(Jv)

        if (present(work)) then
            if (.not.work%validate_size(n, m)) call work%resize(n, m)
            xd => work%xd(1:n)
            fd => work%fd(1:m)
        else
            allocate(x_buffer(n), f_buffer(m))
            xd => x_buffer
            fd => f_buffer
        end if

        do concurrent (i = 1:n)
            xd(i)%x = x(i)
            xd(i)%dx = v(i)
        end do

        call fcn(xd, fd, args = args)

        do i = 1,m
            Jv(i) = fd(i)%dx
        end do
    end subroutine

end module