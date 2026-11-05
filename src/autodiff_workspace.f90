module autodiff_workspace
    !! Reusable storage for dual-valued inputs and outputs.
    use iso_fortran_env
    use autodiff_dual_numbers
    implicit none

    type :: ad_work
        !! Allocatable temporary arrays used during automatic differentiation.
        type(dual), allocatable :: xd(:)
        !! Dual-valued input buffer.
        type(dual), allocatable :: fd(:)
        !! Dual-valued function output buffer.
    end type

contains

    subroutine resize(work, n, m)
        !! Allocate or resize both workspace buffers to the requested dimensions.
        type(ad_work), intent(inout) :: work
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