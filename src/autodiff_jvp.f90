module autodiff_jvp
    !! Forward-mode evaluation of Jacobian-vector products.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces

contains

    subroutine jvp(fcn, x, v, Jv, args)
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
        class(*), intent(inout), optional :: args
            !! Optional caller-owned context passed to the procedure, which 
            !! may update it.

        type(dual), allocatable, dimension(:) :: xd, fd
        integer :: i, n, m

        n = size(x)
        m = size(Jv)

        allocate(xd(n))
        allocate(fd(m))

        do i = 1,n
            xd(i)%x = x(i)
            xd(i)%dx = v(i)
        end do

        call fcn(xd, fd, args = args)

        do i = 1,m
            Jv(i) = fd(i)%dx
        end do
    end subroutine

end module