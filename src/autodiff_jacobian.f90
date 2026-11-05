module autodiff_jacobian
    !! Forward-mode evaluation of vector-function Jacobian matrices.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    implicit none

contains

    subroutine jacobian(fcn, x, J, args)
        !! Compute a vector function's Jacobian by seeding each input basis 
        !! direction.
        procedure(vector_function), pointer, intent(in) :: fcn
            !! User procedure that evaluates the vector-valued function.
        real(real64), intent(in), dimension(:) :: x
            !! Input point at which the Jacobian is evaluated.
        real(real64), intent(out), dimension(:,:) :: J
            !! Jacobian output with shape (number of outputs, number of inputs);
            !! \(J_{ik} = \partial f_i/\partial x_k\).
        class(*), intent(inout), optional :: args
            !! Optional caller-owned context passed to the procedure, which 
            !! may update it.

        type(dual), allocatable, dimension(:) :: xd, fd
        integer :: i, k, m, n

        n = size(x)
        m = size(J,1)

        allocate(xd(n))
        allocate(fd(m))

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