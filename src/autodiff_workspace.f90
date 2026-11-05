module autodiff_workspace
    !! Reusable storage for dual-valued inputs and outputs.
    use iso_fortran_env
    use autodiff_dual_numbers
    implicit none
    private
    public :: workspace

    type :: workspace
        !! Allocatable temporary arrays used during automatic differentiation.
        type(dual), allocatable :: xd(:)
        !! Dual-valued input buffer.
        type(dual), allocatable :: fd(:)
        !! Dual-valued function output buffer.
    contains
        procedure, public :: resize => wrk_resize
            !! Allocates or resizes the workspace buffers.
        procedure, public :: validate_size => wrk_validate_size
            !! Validates the size of the workspace.
    end type

contains

    pure function wrk_validate_size(this, xd_size, fd_size) result(rst)
        !! Validates the size of the workspace arrays is at least
        !! the size requested, if not larger.
        class(workspace), intent(in) :: this
            !! The [[workspace]] object.
        integer, intent(in) :: xd_size
            !! The desired minimum size of 'xd'.
        integer, intent(in) :: fd_size
            !! The desired minimum size of 'fd'.
        logical :: rst
            !! True if both 'xd' and 'fd' meet the minimum required
            !! size; else, false if either of them do not satisfy the
            !! requirements.
        if (.not.allocated(this%xd) .or. .not.allocated(this%fd)) then
            rst = .false.
            return
        end if
        rst = size(this%xd) >= xd_size .and. size(this%fd) >= fd_size
    end function

    subroutine wrk_resize(work, n, m)
        !! Allocate or resize both workspace buffers to the requested dimensions.
        class(workspace), intent(inout) :: work
            !! Workspace whose input and output buffers are managed.
        integer, intent(in) :: n
            !! Required number of input elements.
        integer, intent(in) :: m
            !! Required number of output elements.

        if (.not.allocated(work%xd)) then
            allocate(work%xd(n))
        else if (size(work%xd) /= n) then
            deallocate(work%xd)
            allocate(work%xd(n))
        end if

        if (.not.allocated(work%fd)) then
            allocate(work%fd(m))
        else if (size(work%fd) /= m) then
            deallocate(work%fd)
            allocate(work%fd(m))
        end if
    end subroutine

end module