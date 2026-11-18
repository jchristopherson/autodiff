module autodiff_hyperdual_workspace
    !! Reusable storage for hyper-dual-valued inputs.
    use iso_fortran_env
    use autodiff_hyperdual_numbers
    implicit none
    private
    public :: hyperdual_workspace

    type :: hyperdual_workspace
        !! Allocatable temporary arrays used during second-order automatic
        !! differentiation.
        type(hyperdual), allocatable, dimension(:,:) :: xd_pool
            !! Per-thread hyper-dual-valued input buffers, one column per
            !! thread.  Serial evaluations use the first column.
    contains
        procedure, public :: resize => hwrk_resize
            !! Allocates or grows the workspace buffers.
        procedure, public :: validate_size => hwrk_validate_size
            !! Validates the size of the workspace.
    end type

contains

    function hwrk_validate_size(this, xd_size, nthreads) result(rst)
        !! Validates that the workspace arrays are allocated and at least
        !! the size requested.
        class(hyperdual_workspace), intent(in) :: this
            !! The [[hyperdual_workspace]] object.
        integer, intent(in) :: xd_size
            !! The desired minimum number of input elements per thread.
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
        if (.not.allocated(this%xd_pool)) return
        rst = size(this%xd_pool, 1) >= xd_size .and. &
            size(this%xd_pool, 2) >= nt
    end function

    subroutine hwrk_resize(work, n, nthreads)
        !! Allocates or grows the workspace buffers to at least the requested
        !! dimensions.  Buffers are never shrunk.  When growth is required, a
        !! dimension is sized to the larger of the request and twice its
        !! previous size.  Buffer contents are not preserved on growth.
        class(hyperdual_workspace), intent(inout) :: work
            !! Workspace whose input buffers are managed.
        integer, intent(in) :: n
            !! Required number of input elements per thread.
        integer, intent(in), optional :: nthreads
            !! Required number of per-thread buffer columns.  If omitted,
            !! one column is required.
        integer :: nt, rr, cc

        nt = 1
        if (present(nthreads)) nt = max(1, nthreads)

        if (allocated(work%xd_pool)) then
            if (size(work%xd_pool, 1) >= n .and. &
                size(work%xd_pool, 2) >= nt) return
            rr = size(work%xd_pool, 1)
            cc = size(work%xd_pool, 2)
            if (rr < n) rr = max(n, 2 * rr)
            if (cc < nt) cc = max(nt, 2 * cc)
            deallocate(work%xd_pool)
        else
            rr = n
            cc = nt
        end if
        allocate(work%xd_pool(rr, cc))
    end subroutine

end module
