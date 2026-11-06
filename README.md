# autodiff

A Fortran automatic differentiation library using forward-mode dual numbers.

## Status
[![CMake](https://github.com/jchristopherson/autodiff/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/autodiff/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/autodiff/workflows/fpm/badge.svg)](https://github.com/jchristopherson/autodiff/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/autodiff/).

## Capabilities

| Capability | Description |
| --- | --- |
| Scalar gradients | Compute the gradient of a scalar-valued function with `gradient`. |
| Jacobian matrices | Compute the full Jacobian of a vector-valued function with `jacobian`. |
| Jacobian-vector products | Evaluate `J(x) v` directly with `jvp`, without forming the Jacobian. |
| Dual-number operations | Propagate derivatives through overloaded arithmetic and common intrinsic functions, including trigonometric, exponential, logarithmic, hyperbolic, and error functions. |
| Reusable workspace | Pass a `workspace` to differentiation routines to reuse temporary storage across calls. |
| Optional OpenMP | Evaluate Jacobian columns in parallel when built with OpenMP; control the thread count with `nthreads`. The callback must be thread-safe when using multiple threads. |
| Callback context | Pass optional caller-owned `args` through to the user function or subroutine. |

## Examples

The examples directory contains complete programs for the three main differentiation routines.

### Scalar gradient

[`autodiff_gradient_example.f90`](examples/autodiff_gradient_example.f90) computes the gradient of the Rosenbrock function and compares it with its analytical gradient:

```fortran
call gradient(rosenbrock, x, g)

function rosenbrock(u, args) result(f)
  type(dual), intent(in) :: u(:)
  class(*), intent(in), optional :: args
  type(dual) :: f
  f = (1.0d0 - u(1))**2 + 1.0d2 * (u(2) - u(1)**2)**2
end function
```

The callback operates on `dual` inputs and returns a `dual` result. `gradient` seeds one input direction at a time and returns one partial derivative per element of `x`.

### Vector-function Jacobian

[`autodiff_jacobian_example.f90`](examples/autodiff_jacobian_example.f90) differentiates a two-output system and compares the resulting matrix with the analytical Jacobian:

```fortran
call jacobian(vec_fcn, x, J)

subroutine vec_fcn(u, f, args)
  type(dual), intent(in) :: u(:)
  type(dual), intent(out) :: f(:)
  class(*), intent(in), optional :: args
  f(1) = u(1)**2 + u(2)**2 - 1.0d0
  f(2) = u(1) - u(2)
end subroutine
```

The result has one row per function output and one column per input variable. When the library is built with OpenMP, pass `nthreads` to evaluate columns in parallel.

### Jacobian-vector product

[`autodiff_jvp_example.f90`](examples/autodiff_jvp_example.f90) evaluates the directional derivative of a two-variable model, then checks it against a formed Jacobian multiplied by the same direction:

```fortran
call jvp(vec_fcn, x, v, Jv)

! Inside vec_fcn:
f(1) = a * u(1) - b * u(1) * u(2)
f(2) = c * u(1) * u(2) - d * u(2)
```

`jvp` returns `J(x) v` without allocating or computing the full Jacobian. This is useful when $ J v $ is needed, and it may be expensive or unnecessary to evaluate the full Jacobian and then perform the multiplication.

## Building with CMake

Configure and build the library with a Fortran compiler:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Install it to the selected prefix (the default is platform-dependent):

```sh
cmake --install build --prefix /path/to/install
```

The build defaults to `Release` when no build type is given. Two options
control optional features:

- `AUTODIFF_ENABLE_OPENMP` (default `ON`): links OpenMP when found, so
  `jacobian` can evaluate columns in parallel via its `nthreads` argument.
  Pass `-DAUTODIFF_ENABLE_OPENMP=OFF` to build serially.
- `AUTODIFF_ENABLE_IPO` (default `ON`): enables interprocedural optimization
  (LTO) when the compiler supports it.

The install includes the library, Fortran module files, and CMake package
configuration. A CMake consumer can then use:

```cmake
find_package(autodiff CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE autodiff::autodiff)
```

Include `autodiff` in Fortran source to use the `dual` type
and its overloaded operators and intrinsic functions.

## Building with fpm

Build the library with [fpm](https://github.com/fortran-lang/fpm):

```sh
fpm build
fpm test
```
