module autodiff_jacobian
    !! Forward-mode evaluation of vector-function Jacobian matrices.
    use iso_fortran_env
    use autodiff_dual_numbers
    use autodiff_interfaces
    use autodiff_workspace
    implicit none

contains

    subroutine jacobian(fcn, x, J, work, args, nthreads)
        !! Compute a vector function's Jacobian by seeding each input basis 
        !! direction.
        !!
        !! When 'nthreads' is greater than one and the library is built with
        !! OpenMP, the Jacobian columns are evaluated in parallel.  In this
        !! case, the user procedure must be thread-safe and must not modify
        !! shared state through 'args' without synchronization.
!$ use omp_lib, only : omp_get_thread_num
        procedure(vector_function) :: fcn
            !! User procedure that evaluates the vector-valued function.
        real(real64), intent(in), dimension(:) :: x
            !! Input point at which the Jacobian is evaluated.
        real(real64), intent(out), dimension(:,:) :: J
            !! Jacobian output with shape (number of outputs, number of inputs);
            !! \(J_{ik} = \partial f_i/\partial x_k\).
        class(workspace), intent(inout), target, optional :: work
            !! An optional workspace to use in order to avoid internal memory
            !! allocations.  If not sized properly, this routine will grow
            !! the workspace to suit its needs.
        class(*), intent(in), optional :: args
            !! Optional caller-owned context passed to the procedure.
        integer, intent(in), optional :: nthreads
            !! The maximum number of threads to use when evaluating the
            !! Jacobian columns.  The default is 1 (serial).

        type(dual), allocatable, target, dimension(:,:) :: x_buffer, f_buffer
        type(dual), pointer, dimension(:,:) :: xd, fd
        integer :: i, k, m, n, nt, tid

        n = size(x)
        m = size(J,1)
        nt = 1
        if (present(nthreads)) nt = max(1, nthreads)

        if (present(work)) then
            if (.not.work%validate_size(n, m, nt)) call work%resize(n, m, nt)
            xd => work%xd_pool(1:n, 1:nt)
            fd => work%fd_pool(1:m, 1:nt)
        else
            allocate(x_buffer(n, nt), f_buffer(m, nt))
            xd => x_buffer
            fd => f_buffer
        end if

!$omp parallel num_threads(nt) default(shared) private(i, k, tid)
        tid = 1
!$      tid = omp_get_thread_num() + 1

        ! Each thread owns one column of the buffers.  Inputs are set once
        ! per thread, and only the seeded direction is toggled per column.
        do i = 1, n
            xd(i,tid)%x = x(i)
            xd(i,tid)%dx = 0.0d0
        end do

!$omp do schedule(static)
        do k = 1, n
            xd(k,tid)%dx = 1.0d0

            call fcn(xd(:,tid), fd(:,tid), args = args)

            do i = 1, m
                J(i,k) = fd(i,tid)%dx
            end do

            xd(k,tid)%dx = 0.0d0
        end do
!$omp end do
!$omp end parallel

    end subroutine


end module

! module autodiff_jacobian
!     !! Forward-mode evaluation of vector-function Jacobian matrices.
!     use iso_fortran_env
!     use autodiff_dual_numbers
!     use autodiff_interfaces
!     use autodiff_workspace
!     implicit none
!     private
!     public :: jacobian

! contains

!     subroutine jacobian(fcn, x, J, work, args, nthreads)
!         !! Compute a vector function's Jacobian by seeding each input basis 
!         !! direction.
!         !!
!         !! When 'nthreads' is greater than one and the library is built with
!         !! OpenMP, the Jacobian columns are evaluated in parallel.  In this
!         !! case, the user procedure must be thread-safe and must not modify
!         !! shared state through 'args' without synchronization.
!         procedure(vector_function) :: fcn
!             !! User procedure that evaluates the vector-valued function.
!         real(real64), intent(in), dimension(:) :: x
!             !! Input point at which the Jacobian is evaluated.
!         real(real64), intent(out), dimension(:,:) :: J
!             !! Jacobian output with shape (number of outputs, number of inputs);
!             !! \(J_{ik} = \partial f_i/\partial x_k\).
!         class(workspace), intent(inout), target, optional :: work
!             !! An optional workspace to use in order to avoid internal memory
!             !! allocations.  If not sized properly, this routine will grow
!             !! the workspace to suit its needs.
!         class(*), intent(in), optional :: args
!             !! Optional caller-owned context passed to the procedure.
!         integer, intent(in), optional :: nthreads
!             !! The maximum number of threads to use when evaluating the
!             !! Jacobian columns.  The default is 1 (serial).

!         type(dual), allocatable, dimension(:,:) :: x_buffer, f_buffer
!         integer :: n, m, nt

!         n = size(x)
!         m = size(J,1)
!         nt = 1
!         if (present(nthreads)) nt = max(1, nthreads)

!         if (present(work)) then
!             if (.not.work%validate_size(n, m, nt)) call work%resize(n, m, nt)
!             if (nt > 1) then
!                 call jacobian_parallel(fcn, x, J, work%xd_pool(:, 1:nt), &
!                     work%fd_pool(:, 1:nt), args = args)
!             else
!                 call jacobian_serial(fcn, x, J, work%xd_pool(:, 1:nt), &
!                     work%fd_pool(:, 1:nt), args = args)
!             end if
!         else
!             allocate(x_buffer(n, nt), f_buffer(m, nt))
!             if (nt > 1) then
!                 call jacobian_parallel(fcn, x, J, x_buffer, f_buffer, args = args)
!             else
!                 call jacobian_serial(fcn, x, J, x_buffer, f_buffer, args = args)
!             end if
!         end if

!     end subroutine

!     subroutine jacobian_serial(fcn, x, J, xd, fd, args)
!         !! Serial Jacobian evaluation using the first column of the buffers.
!         procedure(vector_function) :: fcn
!             !! User procedure that evaluates the vector-valued function.
!         real(real64), intent(in), dimension(:) :: x
!             !! Input point at which the Jacobian is evaluated.
!         real(real64), intent(out), dimension(:,:) :: J
!             !! Jacobian output with shape (number of outputs, number of inputs).
!         type(dual), intent(inout), dimension(:,:) :: xd, fd
!             !! Input and output buffers; only the first column is used.
!         class(*), intent(in), optional :: args
!             !! Optional caller-owned context passed to the procedure.

!         integer :: i, k, m, n

!         n = size(x)
!         m = size(J,1)

!         do i = 1, n
!             xd(i,1)%x = x(i)
!             xd(i,1)%dx = 0.0d0
!         end do

!         do k = 1, n
!             xd(k,1)%dx = 1.0d0

!             if (present(args)) then
!                 call fcn(xd(1:n,1), fd(1:m,1), args = args)
!             else
!                 call fcn(xd(1:n,1), fd(1:m,1))
!             end if

!             do i = 1, m
!                 J(i,k) = fd(i,1)%dx
!             end do

!             xd(k,1)%dx = 0.0d0
!         end do

!     end subroutine

!     subroutine jacobian_parallel(fcn, x, J, xd, fd, args)
!         !! OpenMP Jacobian evaluation; each thread owns one column of the buffers.
! !$ use omp_lib, only : omp_get_thread_num
!         procedure(vector_function) :: fcn
!             !! User procedure that evaluates the vector-valued function.
!         real(real64), intent(in), dimension(:) :: x
!             !! Input point at which the Jacobian is evaluated.
!         real(real64), intent(out), dimension(:,:) :: J
!             !! Jacobian output with shape (number of outputs, number of inputs).
!         type(dual), intent(inout), dimension(:,:) :: xd, fd
!             !! Per-thread input and output buffers, one column per thread.
!         class(*), intent(in), optional :: args
!             !! Optional caller-owned context passed to the procedure.

!         integer :: i, k, m, n, nt, tid

!         n = size(x)
!         m = size(J,1)
!         nt = size(xd, 2)

! !$omp parallel num_threads(nt) default(shared) private(i, k, tid)
!         tid = 1
! !$      tid = omp_get_thread_num() + 1

!         ! Each thread owns one column of the buffers.  Inputs are set once
!         ! per thread, and only the seeded direction is toggled per column.
!         do i = 1, n
!             xd(i,tid)%x = x(i)
!             xd(i,tid)%dx = 0.0d0
!         end do

! !$omp do schedule(static)
!         do k = 1, n
!             xd(k,tid)%dx = 1.0d0

!             if (present(args)) then
!                 call fcn(xd(1:n,tid), fd(1:m,tid), args = args)
!             else
!                 call fcn(xd(1:n,tid), fd(1:m,tid))
!             end if

!             do i = 1, m
!                 J(i,k) = fd(i,tid)%dx
!             end do

!             xd(k,tid)%dx = 0.0d0
!         end do
! !$omp end do
! !$omp end parallel

!     end subroutine

! end module