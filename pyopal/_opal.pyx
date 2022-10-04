# distutils: language = c++
# cython: language_level=3, linetrace=True, embedsignature=True, binding=True
"""Bindings to Opal, a SIMD-accelerated pairwise sequence aligner.

References:
    - Korpar M., Šošić M., Blažeka D., Šikić M. 
      *SW#db: GPU-Accelerated Exact Sequence Similarity Database Search*. 
      PLoS One. 2015;10(12):e0145857. Published 2015 Dec 31. 
      :doi:`10.1371/journal.pone.0145857`.

"""

# --- C imports ----------------------------------------------------------------

cimport opal


# --- Runtime CPU detection ----------------------------------------------------

_TARGET_CPU           = TARGET_CPU
_SSE2_RUNTIME_SUPPORT = False
_SSE2_BUILD_SUPPORT   = False
_SSE4_RUNTIME_SUPPORT = False
_SSE4_BUILD_SUPPORT   = False
_AVX2_RUNTIME_SUPPORT = False
_AVX2_BUILD_SUPPORT   = False
_NEON_RUNTIME_SUPPORT = False
_NEON_BUILD_SUPPORT   = False

IF TARGET_CPU == "x86" and TARGET_SYSTEM in ("freebsd", "linux_or_android", "macos", "windows"):
    from cpu_features.x86 cimport GetX86Info, X86Info
    cdef X86Info cpu_info = GetX86Info()
    _SSE2_BUILD_SUPPORT   = SSE2_BUILD_SUPPORT
    _SSE2_RUNTIME_SUPPORT = SSE2_BUILD_SUPPORT and cpu_info.features.sse2 != 0
    _SSE4_BUILD_SUPPORT   = SSE4_BUILD_SUPPORT
    _SSE4_RUNTIME_SUPPORT = SSE4_BUILD_SUPPORT and cpu_info.features.sse4_1 != 0
    _AVX2_BUILD_SUPPORT   = AVX2_BUILD_SUPPORT
    _AVX2_RUNTIME_SUPPORT = AVX2_BUILD_SUPPORT and cpu_info.features.avx2 != 0
ELIF TARGET_CPU == "arm":
    from cpu_features.arm cimport GetArmInfo, ArmInfo
    cdef ArmInfo arm_info = GetArmInfo()
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT and arm_info.features.neon != 0
ELIF TARGET_CPU == "aarch64":
    _NEON_BUILD_SUPPORT   = NEON_BUILD_SUPPORT
    _NEON_RUNTIME_SUPPORT = NEON_BUILD_SUPPORT # always runtime support on Aarch64


if _SSE2_RUNTIME_SUPPORT:
    _BEST_BACKEND = simd_backend.SSE2
elif _NEON_RUNTIME_SUPPORT:
    _BEST_BACKEND = simd_backend.NEON
else:
    _BEST_BACKEND = simd_backend.GENERIC

cdef enum simd_backend:
    NONE = 0
    GENERIC = 1
    SSE2 = 2
    NEON = 3



IF SSE2_BUILD_SUPPORT:
    from pyopal._opal_sse2 cimport opalSearchDatabaseSSE2
IF AVX2_BUILD_SUPPORT:
    from pyopal._opal_avx2 cimport opalSearchDatabaseAVX2


cdef class SearchResults:
    pass

cdef class SearchResult:
    pass


cdef class Database:

    def __init__(self, *sequences):
        pass
    
    def search(self, bytes sequence):
        if _AVX2_RUNTIME_SUPPORT:
            opalSearchDatabaseAVX2(
                NULL,
                0,
                NULL,
                0,
                NULL,
                0,
                0,
                NULL,
                0,
                NULL,
                0,
                0,
                0,
            )
        elif _SSE2_RUNTIME_SUPPORT:
            opalSearchDatabaseSSE2(
                NULL,
                0,
                NULL,
                0,
                NULL,
                0,
                0,
                NULL,
                0,
                NULL,
                0,
                0,
                0,
            )
    