set(CMAKE_C_COMPILER "/opt/ecrc/gcc/7.2.0/ub16/bin/gcc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "GNU")
set(CMAKE_C_COMPILER_VERSION "7.2.0")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "")
set(CMAKE_C_SIMULATE_VERSION "")



set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "/opt/ecrc/gcc/7.2.0/ub16/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "/opt/ecrc/gcc/7.2.0/ub16/bin/gcc-ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_COMPILER_IS_GNUCC 1)
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "x86_64-linux-gnu")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "gcc;gcc_s;c;gcc;gcc_s")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/opt/ecrc/gcc/7.2.0/ub16/lib64;/opt/ecrc/gcc/7.2.0/ub16/lib/gcc/x86_64-linux-gnu/7.2.0;/lib/x86_64-linux-gnu;/lib64;/usr/lib/x86_64-linux-gnu;/opt/ecrc/r-base/3.5.1-mkl/ub16/lib/R/lib;/opt/ecrc/r-base/3.5.1-mkl/ub16/lib;/opt/ecrc/pcre/8.40/ub16/lib;/opt/ecrc/starpu/1.2.4-gcc-7.2.0-mkl-openmpi-3.0.0-fxt/ub16/lib;/opt/ecrc/fxt/0.3.7-gcc-7.2.0/ub16/lib;/opt/ecrc/hwloc/1.11.8-gcc-7.2.0/ub16/lib;/opt/ecrc/nlopt/2.4.2-gcc-7.2.0/ub16/lib;/opt/ecrc/gsl/2.4-gcc-7.2.0/ub16/lib;/opt/ecrc/mkl/2018-update-1/mkl/lib/intel64;/opt/ecrc/mkl/2018-update-1/lib/intel64;/opt/ecrc/gcc/7.2.0/ub16/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
