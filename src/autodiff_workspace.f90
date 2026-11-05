module autodiff_workspace
    !! Reusable storage for dual-valued inputs and outputs.
    use iso_fortran_env
    use autodiff_dual_numbers
    implicit none
    private
    public :: workspace

    type :: workspace
        !! Allocatable temporary arrays used during automatic differentiation.
        type(dual), allocatable, dimension(:) :: xd
        !! Dual-valued input buffer.
        type(dual), allocatable, dimension(:) :: fd
        !! Dual-valued function output buffer.
        type(dual), allocatable, dimension(:,:) :: xd_pool
        !! Per-thread dual-valued input buffers, one column per thread.
        type(dual), allocatable, dimension(:,:) :: fd_pool
        !! Per-thread dual-valued function output buffers, one column per thread.
    contains
        procedure, public :: resize => wrk_resize
            !! Allocates or grows the workspace buffers.
        procedure, public :: validate_size => wrk_validate_size
            !! Validates the size of the workspace.
    end type

contains

    function wrk_validate_size(this, xd_size, fd_size, nthreads) result(rst)
        !! Validates that the workspace arrays are allocated and at least
        !! the size requested.
        class(workspace), intent(in) :: this
            !! The [[workspace]] object.
        integer, intent(in) :: xd_size
            !! The desired minimum size of 'xd'.
        integer, intent(in) :: fd_size
            !! The desired minimum size of 'fd'.
        integer, intent(in), optional :: nthreads
            !! The desired number of per-thread buffer columns.  If omitted,
            !! one column is required.
        logical :: rst
            !! True if all buffers meet the minimum required size; else,
            !! false.
        integer :: nt

        nt = 1
        if (present(nthreads)) nt = max(1, nthreads)

        rst = .false.
        if (.not.allocated(this%xd)) return
        if (.not.allocated(this%fd)) return
        if (.not.allocated(this%xd_pool)) return
        if (.not.allocated(this%fd_pool)) return

        rst = size(this%xd) >= xd_size .and. size(this%fd) >= fd_size
        if (.not.rst) return
        rst = size(this%xd_pool, 1) >= xd_size .and. &
            size(this%xd_pool, 2) >= nt .and. &
            size(this%fd_pool, 1) >= fd_size .and. &
            size(this%fd_pool, 2) >= nt
    end function

    subroutine wrk_resize(work, n, m, nthreads)
        !! Allocates or grows the workspace buffers to at least the requested
        !! dimensions.  Buffers are never shrunk.  When growth is required, a
        !! buffer is sized to the larger of the request and twice its
        !! previous size.  Buffer contents are not preserved on growth.
        class(workspace), intent(inout) :: work
            !! Workspace whose input and output buffers are managed.
        integer, intent(in) :: n
            !! Required number of input elements.
        integer, intent(in) :: m
            !! Required number of output elements.
        integer, intent(in), optional :: nthreads
            !! Required number of per-thread buffer columns.  If omitted,
            !! one column is required.
        integer :: nt

        nt = 1
        if (present(nthreads)) nt = max(1, nthreads)

        call grow_vector(work%xd, n)
        call grow_vector(work%fd, m)
        call grow_matrix(work%xd_pool, n, nt)
        call grow_matrix(work%fd_pool, m, nt)
    end subroutine

    subroutine grow_vector(v, n)
        !! Grows an allocatable vector to at least 'n' elements.
        type(dual), allocatable, intent(inout) :: v(:)
            !! The vector to grow.
        integer, intent(in) :: n
            !! The required number of elements.
        integer :: nn

        if (allocated(v)) then
            if (size(v) >= n) return
            nn = max(n, 2 * size(v))
            deallocate(v)
        else
            nn = n
        end if
        allocate(v(nn))
    end subroutine

    subroutine grow_matrix(a, nr, nc)
        !! Grows an allocatable matrix to at least 'nr' rows and 'nc' columns.
        type(dual), allocatable, intent(inout) :: a(:,:)
            !! The matrix to grow.
        integer, intent(in) :: nr
            !! The required number of rows.
        integer, intent(in) :: nc
            !! The required number of columns.
        integer :: rr, cc

        if (allocated(a)) then
            if (size(a, 1) >= nr .and. size(a, 2) >= nc) return
            rr = size(a, 1)
            cc = size(a, 2)
            if (rr < nr) rr = max(nr, 2 * rr)
            if (cc < nc) cc = max(nc, 2 * cc)
            deallocate(a)
        else
            rr = nr
            cc = nc
        end if
        allocate(a(rr, cc))
    end subroutine

end module