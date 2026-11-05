# autodiff

A Fortran automatic differentiation library.

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

The install includes the library, Fortran module files, and CMake package
configuration. A CMake consumer can then use:

```cmake
find_package(autodiff CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE autodiff::autodiff)
```

Include `autodiff_dual_numbers` in Fortran source to use the `real_ad` type
and its overloaded operators and intrinsic functions.

## Building with fpm

Build the library with [fpm](https://github.com/fortran-lang/fpm):

```sh
fpm build
fpm test
```
