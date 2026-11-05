module autodiff_interfaces
    !! Abstract procedure interfaces used by the forward-mode differentiation routines.
    use iso_fortran_env
    use autodiff_dual_numbers
    implicit none
    private
    public :: scalar_function
    public :: vector_function

    abstract interface

        function scalar_function(x, args) result(f)
            !! Evaluate a scalar function on dual inputs; the result's `dx` is 
            !! the seeded directional derivative.
            use autodiff_dual_numbers, only : dual
                !! Provides the dual-number argument and result type.
            type(dual), intent(in), dimension(:) :: x
                !! Input values and derivative seeds; for seed vector \(v\), the
                !! returned derivative is \(D f(x)[v] = \nabla f(x)^{T}v\).
            class(*), intent(in), optional :: args
                !! Optional caller-owned context passed to the function.  It is
                !! read-only from the function's perspective.
            type(dual) :: f
                !! Function value in `f%x` and directional derivative in `f%dx`.
        end function

        subroutine vector_function(x, f, args)
            !! Evaluate a vector function on dual inputs; each output's `dx` is 
            !! its directional derivative.
            use autodiff_dual_numbers, only : dual
                !! Provides the dual-number input and output type.
            type(dual), intent(in), dimension(:) :: x
                !! Input values and derivative seeds; for seed vector \(v\), 
                !! outputs propagate the direction \(v\).
            type(dual), intent(out), dimension(:) :: f
                !! Output values and derivatives; the derivatives form 
                !! \(J_f(x)v\), where \(J_f(x)\) is the Jacobian.
            class(*), intent(in), optional :: args
                !! Optional caller-owned context passed to the procedure.  It is
                !! read-only from the procedure's perspective.
        end subroutine

    end interface

end module