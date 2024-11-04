get_property(OPAL_CPP GLOBAL PROPERTY OPAL_CPP)

list(APPEND _PLATFORMS avx2 sse2 sse4 neon)
file(READ pyx.in TEMPLATE)

set(AVX2_DEFINE_MACROS "-D__AVX2__")
set(SSE2_DEFINE_MACROS "-D__SSE2__")
set(SSE4_DEFINE_MACROS "-D__SSE4_1__")
set(NEON_DEFINE_MACROS "-D__ARM_NEON__")

foreach(_simd IN LISTS _PLATFORMS)
    string(TOUPPER ${_simd} SIMD)
    if(HAVE_${SIMD})
        # copy `opal.cpp` to a new file to force CMake to rebuild it with
        # different compile flags (SIMD-specific)
        add_custom_command(
            OUTPUT
                opal_${_simd}.cpp
            COMMAND
                cmake -E copy ${PROJECT_SOURCE_DIR}/vendor/opal/src/opal.cpp opal_${_simd}.cpp
            DEPENDS
                ${PROJECT_SOURCE_DIR}/vendor/opal/src/opal.cpp
        )
        # generate the Cython extension code from the template
        file(CONFIGURE CONTENT ${TEMPLATE} OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/${_simd}.pyx)
        cython_extension(${_simd} LINKS opal EXTRA_SOURCES opal_${_simd}.cpp)
        # add SIMD flags to the extension
        string(REPLACE " " ";" IMPL_FLAGS ${${SIMD}_C_FLAGS})
        foreach(_flag IN LISTS IMPL_FLAGS)
           target_compile_options(pyopal.platform.${_simd} PRIVATE ${_flag})
           target_compile_definitions(pyopal.platform.${_simd} PRIVATE ${${SIMD}_DEFINE_MACROS})
        endforeach()
    endif()
endforeach()