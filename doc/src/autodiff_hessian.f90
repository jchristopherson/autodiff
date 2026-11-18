module autodiff_hessian
    !! Second-order forward-mode evaluation of scalar-function gradients and Hessians.
    use iso_fortran_env
    use autodiff_hyperdual_numbers
    use autodiff_interfaces
    use autodiff_hyperdual_workspace
    implicit none
    private
    public :: hessian
    public :: hvp

contains
! ------------------------------------------------------------------------------
    subroutine hessian(fcn, x, H, g, f, work, args, nthreads)
        !! Compute a scalar function's Hessian, and optionally its gradient and
        !! value, by seeding hyper-dual inputs.  The Hessian is symmetric and
        !! requires \(n(n+1)/2\) evaluations of `fcn`; the gradient and value
        !! are taken from the diagonal evaluations at no additional cost.
        !!
        !! When 'nthreads' is greater than one and the library is built with
        !! OpenMP, the Hessian rows are evaluated in parallel.  In this case,
        !! the user function must be thread-safe and must not modify shared
        !! state through 'args' without synchronization.
!$ use omp_lib, only : omp_get_thread_num
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
        class(hyperdual_workspace), intent(inout), target, optional :: work
            !! An optional workspace to use in order to avoid internal memory
            !! allocations.  If not sized properly, this routine will grow
            !! the workspace to suit its needs.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the function.
        integer, intent(in), optional :: nthreads
            !! The maximum number of threads to use when evaluating the
            !! Hessian rows.  The default is 1 (serial).

        type(hyperdual), allocatable, target, dimension(:,:) :: x_buffer
        type(hyperdual), pointer, dimension(:,:) :: xd
        type(hyperdual) :: fd
        integer :: i, n, nt, tid

        n = size(x)
        nt = 1
        if (present(nthreads)) nt = max(1, nthreads)

        if (present(work)) then
            if (.not.work%validate_size(n, nt)) call work%resize(n, nt)
            xd => work%xd_pool(1:n, 1:nt)
        else
            allocate(x_buffer(n, nt))
            xd => x_buffer
        end if

        if (n == 0) then
            if (present(f)) then
                fd = fcn(xd(:,1), args = args)
                f = fd%x
            end if
            return
        end if

        if (nt == 1) then
            ! Run serial, but do not engage with OpenMP
            call init_inputs(x, xd(:,1))
            do i = 1, n
                call hessian_row(fcn, xd(:,1), i, H, g, f, args)
            end do
        else
            ! Utilize OpenMP.  Row i requires n - i + 1 evaluations, so rows
            ! are dealt cyclically to balance the triangular workload.
        !$omp parallel num_threads(nt) default(shared) private(i, tid)
            tid = 1
        !$      tid = omp_get_thread_num() + 1

            ! Each thread owns one column of the buffer.  Inputs are set once
            ! per thread, and only the seeded directions are toggled per entry.
            call init_inputs(x, xd(:,tid))

        !$omp do schedule(static, 1)
            do i = 1, n
                call hessian_row(fcn, xd(:,tid), i, H, g, f, args)
            end do
        !$omp end do
        !$omp end parallel
        end if
    end subroutine

! ------------------------------------------------------------------------------
    subroutine hvp(fcn, x, v, Hv, work, args)
        !! Evaluate the product of a scalar function's Hessian and a vector
        !! without forming the Hessian.  Each evaluation seeds the first
        !! derivative direction with `v` and the second with a unit vector,
        !! so the mixed derivative yields one component of \(Hv\).  This
        !! requires `size(x)` evaluations of `fcn`.
        procedure(hyperdual_function) :: fcn
            !! User function whose result carries the value and derivatives.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which the Hessian is evaluated.
        real(real64), intent(in), dimension(:) :: v
            !! Direction vector with one component per input.
        real(real64), intent(out), dimension(:) :: Hv
            !! Product of the Hessian and v, satisfying \(Hv = H_f(x)v\).
        class(hyperdual_workspace), intent(inout), target, optional :: work
            !! An optional workspace to use in order to avoid internal memory
            !! allocations.  If not sized properly, this routine will grow
            !! the workspace to suit its needs.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the function.

        type(hyperdual), allocatable, target, dimension(:,:) :: x_buffer
        type(hyperdual), pointer, dimension(:,:) :: xd
        type(hyperdual) :: fd
        integer :: i, n

        n = size(x)
        if (present(work)) then
            if (.not.work%validate_size(n)) call work%resize(n)
            xd => work%xd_pool(1:n, 1:1)
        else
            allocate(x_buffer(n, 1))
            xd => x_buffer
        end if

        call init_inputs(x, xd(:,1))
        do i = 1, n
            xd(i,1)%dx = v(i)
        end do

        do i = 1, n
            xd(i,1)%dy = 1.0d0
            fd = fcn(xd(:,1), args = args)
            Hv(i) = fd%dxy
            xd(i,1)%dy = 0.0d0
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine init_inputs(x, xd)
        !! Copies the evaluation point into a hyper-dual buffer and clears
        !! all derivative seeds.
        real(real64), intent(in), dimension(:) :: x
            !! Point at which derivatives are evaluated.
        type(hyperdual), intent(inout), dimension(:) :: xd
            !! Hyper-dual input buffer with at least size(x) elements.
        integer :: k

        do k = 1, size(x)
            xd(k) = x(k)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    subroutine hessian_row(fcn, xd, i, H, g, f, args)
        !! Evaluates the upper-triangular entries of row 'i' of the Hessian
        !! and mirrors them into column 'i'.  The diagonal evaluation also
        !! supplies the gradient component and, for the first row, the value.
        !! The seeds in 'xd' are restored to zero on return.
        procedure(hyperdual_function) :: fcn
            !! User function whose result carries the value and derivatives.
        type(hyperdual), intent(inout), dimension(:) :: xd
            !! Hyper-dual input buffer with values set and seeds cleared.
        integer, intent(in) :: i
            !! Row index.
        real(real64), intent(inout), dimension(:,:) :: H
            !! Hessian output.
        real(real64), intent(inout), optional, dimension(:) :: g
            !! Optional gradient output.
        real(real64), intent(inout), optional :: f
            !! Optional function value output.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the function.

        type(hyperdual) :: fd
        integer :: j

        xd(i)%dx = 1.0d0
        do j = i, size(xd)
            xd(j)%dy = 1.0d0
            fd = fcn(xd, args = args)
            H(i, j) = fd%dxy
            H(j, i) = fd%dxy
            if (j == i) then
                if (present(g)) g(i) = fd%dx
                if (i == 1 .and. present(f)) f = fd%x
            end if
            xd(j)%dy = 0.0d0
        end do
        xd(i)%dx = 0.0d0
    end subroutine

! ------------------------------------------------------------------------------
end module
