import configparser
import functools
import glob
import itertools
import io
import multiprocessing.pool
import os
import platform
import re
import setuptools
import setuptools.extension
import subprocess
import sys
from distutils.command.clean import clean as _clean
from distutils.errors import CompileError
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist
from setuptools.extension import Library

try:
    from Cython.Build import cythonize
except ImportError as err:
    cythonize = err

# --- Extension with SIMD requirements -----------------------------------------

class Extension(setuptools.extension.Extension):
    def __init__(self, *args, **kwargs):
        self._needs_stub = False
        self.requires = kwargs.pop("requires", None)
        super().__init__(*args, **kwargs)


# --- Constants -----------------------------------------------------------------

SETUP_FOLDER = os.path.realpath(os.path.join(__file__, os.pardir))
INCLUDE_FOLDER = os.path.join(SETUP_FOLDER, "vendor", "opal", "src")

MACHINE = platform.machine()
if re.match("^mips", MACHINE):
    TARGET_CPU = "mips"
elif re.match("^(aarch64|arm64)$", MACHINE):
    TARGET_CPU = "aarch64"
elif re.match("^arm", MACHINE):
    TARGET_CPU = "arm"
elif re.match("(x86_64)|(AMD64|amd64)|(^i.86$)", MACHINE):
    TARGET_CPU = "x86"
elif re.match("^(powerpc|ppc)", MACHINE):
    TARGET_CPU = "ppc"
else:
    TARGET_CPU = None

SYSTEM = platform.system()
if SYSTEM == "Linux" or SYSTEM == "Java":
    TARGET_SYSTEM = "linux_or_android"
elif SYSTEM.endswith("FreeBSD"):
    TARGET_SYSTEM = "freebsd"
elif SYSTEM == "Darwin":
    TARGET_SYSTEM = "macos"
elif SYSTEM.startswith(("Windows", "MSYS", "MINGW", "CYGWIN")):
    TARGET_SYSTEM = "windows"
else:
    TARGET_SYSTEM = None

# --- Utils ------------------------------------------------------------------

_HEADER_PATTERN = re.compile("^@@ -(\d+),?(\d+)? \+(\d+),?(\d+)? @@$")


def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def _patch_osx_compiler(compiler):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of SSE2.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next(
            (
                i
                for i in range(1, len(flags))
                if flags[i - 1] == "-arch" and flags[i] != platform.machine()
            ),
            None,
        )
        if i is not None:
            flags.pop(i)
            flags.pop(i - 1)


# --- Commands ------------------------------------------------------------------


class sdist(_sdist):
    """A `sdist` that generates a `pyproject.toml` on the fly."""

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", "build-backend", '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)
        # run the rest of the packaging
        _sdist.run(self)


class build_ext(_build_ext):
    """A `build_ext` that adds various SIMD flags and defines."""

    # --- Compatibility with `setuptools.Command`

    user_options = _build_ext.user_options + [
        (
            "disable-avx2",
            None,
            "Force compiling the extension without AVX2 instructions",
        ),
        (
            "disable-sse2",
            None,
            "Force compiling the extension without SSE2 instructions",
        ),
        (
            "disable-sse4",
            None,
            "Force compiling the extension without SSE4.1 instructions",
        ),
        (
            "disable-neon",
            None,
            "Force compiling the extension without NEON instructions",
        ),
    ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.disable_avx2 = False
        self.disable_sse2 = False
        self.disable_sse4 = False
        self.disable_neon = False

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # record SIMD-specific options
        self._simd_supported = dict(AVX2=False, SSE2=False, NEON=False, SSE4=False)
        self._simd_defines = dict(AVX2=[], SSE2=[], NEON=[], SSE4=[])
        self._simd_flags = dict(AVX2=[], SSE2=[], NEON=[], SSE4=[])
        self._simd_disabled = {
            "AVX2": self.disable_avx2,
            "SSE2": self.disable_sse2,
            "SSE4": self.disable_sse4,
            "NEON": self.disable_neon,
        }
        # transfer arguments to the build_clib method
        self._clib_cmd = self.get_finalized_command("build_clib")
        self._clib_cmd.debug = self.debug
        self._clib_cmd.force = self.force
        self._clib_cmd.verbose = self.verbose
        self._clib_cmd.define = self.define
        self._clib_cmd.include_dirs = self.include_dirs
        self._clib_cmd.compiler = self.compiler
        self._clib_cmd.parallel = self.parallel

    # --- Autotools-like helpers ---

    def _check_simd_generic(self, name, flags, header, vector, set, op, extract):
        _eprint("checking whether compiler can build", name, "code", end="... ")

        base = "have_{}".format(name)
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        binfile = self.compiler.executable_filename(base, output_dir=self.build_temp)
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write(
                """
                #include <{}>
                int main() {{
                    {}      a = {}(1);
                            a = {}(a);
                    short   x = {}(a, 1);
                    return (x == 1) ? 0 : 1;
                }}
            """.format(
                    header, vector, set, op, extract
                )
            )

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_preargs=flags)
            self.compiler.link_executable(objects, base, extra_preargs=flags, output_dir=self.build_temp)
            subprocess.run([binfile], check=True)
        except CompileError:
            _eprint("no")
            return False
        except subprocess.CalledProcessError:
            _eprint("yes, but cannot run code")
            return True  # assume we are cross-compiling, and still build
        else:
            if not flags:
                _eprint("yes")
            else:
                _eprint("yes, with {}".format(" ".join(flags)))
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)
            if os.path.isfile(binfile):
                os.remove(binfile)

    def _avx2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:AVX2"]
        return ["-mavx", "-mavx2"]

    def _check_avx2(self):
        return self._check_simd_generic(
            "AVX2",
            self._avx2_flags(),
            header="immintrin.h",
            vector="__m256i",
            set="_mm256_set1_epi16",
            op="_mm256_abs_epi32",
            extract="_mm256_extract_epi16",
        )

    def _sse2_flags(self):
        if self.compiler.compiler_type == "msvc":
            return []
        return ["-msse2"]

    def _check_sse2(self):
        return self._check_simd_generic(
            "SSE2",
            self._sse2_flags(),
            header="emmintrin.h",
            vector="__m128i",
            set="_mm_set1_epi16",
            op="_mm_move_epi64",
            extract="_mm_extract_epi16",
        )

    def _sse4_flags(self):
        if self.compiler.compiler_type == "msvc":
            return ["/arch:AVX"]
        return ["-msse4.1"]

    def _check_sse4(self):
        return self._check_simd_generic(
            "SSE4",
            self._sse4_flags(),
            header="smmintrin.h",
            vector="__m128i",
            set="_mm_set1_epi8",
            op="_mm_cvtepi8_epi16",
            extract="_mm_extract_epi16",
        )

    def _neon_flags(self):
        return ["-mfpu=neon"] if TARGET_CPU == "arm" else []

    def _check_neon(self):
        return self._check_simd_generic(
            "NEON",
            self._neon_flags(),
            header="arm_neon.h",
            vector="int16x8_t",
            set="vdupq_n_s16",
            op="vabsq_s16",
            extract="vgetq_lane_s16",
        )

    def _check_getid(self):
        _eprint('checking whether `PyInterpreterState_GetID` is available')

        base = "have_getid"
        testfile = os.path.join(self.build_temp, "{}.c".format(base))
        objects = []

        self.mkpath(self.build_temp)
        with open(testfile, "w") as f:
            f.write("""
            #include <stdint.h>
            #include <stdlib.h>
            #include <Python.h>

            int main(int argc, char *argv[]) {{
                PyInterpreterState_GetID(NULL);
                return 0;
            }}
            """)

        if self.compiler.compiler_type == "msvc":
            flags = ["/WX"]
        else:
            flags = ["-Werror=implicit-function-declaration"]

        try:
            self.mkpath(self.build_temp)
            objects = self.compiler.compile([testfile], extra_postargs=flags)
        except CompileError:
            _eprint("no")
            return False
        else:
            _eprint("yes")
            return True
        finally:
            os.remove(testfile)
            for obj in filter(os.path.isfile, objects):
                os.remove(obj)

    # --- Build code ---

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "with", self.compiler.compiler_type, "compiler")

        # add debug symbols if we are building in debug mode
        if self.debug:
            if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
                ext.extra_compile_args.append("-g")
            elif self.compiler.compiler_type == "msvc":
                ext.extra_compile_args.append("/Z7")
            if sys.implementation.name == "cpython":
                ext.define_macros.append(("CYTHON_TRACE_NOGIL", 1))
        else:
            ext.define_macros.append(("CYTHON_WITHOUT_ASSERTIONS", 1))

        # add C++17 flags (<shared_mutex>) and silence some warnings
        if self.compiler.compiler_type in {"unix", "cygwin", "mingw32"}:
            ext.extra_compile_args.extend([
                "-funroll-loops",
                "-std=c++17",
                "-Wno-unused-variable",
                "-Wno-maybe-uninitialized",
                "-Wno-return-type"
            ])
        elif self.compiler.compiler_type == "msvc":
            ext.extra_compile_args.append("/std:c17")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            ext.define_macros.append(("WIN32", 1))

        # update link and include directories
        for name in ext.libraries:
            lib = self._clib_cmd.get_library(name)
            libfile = self.compiler.library_filename(
                lib.name, output_dir=self._clib_cmd.build_clib
            )
            ext.depends.append(libfile)
            ext.extra_objects.append(libfile)

        # check if `PyInterpreterState_GetID` is defined
        if self._check_getid():
            ext.define_macros.append(("HAS_PYINTERPRETERSTATE_GETID", 1))

        # build the rest of the extension as normal
        ext._needs_stub = False

        # compile extension in its own folder: since we need to compile
        # `opal.cpp` several times with different flags, we cannot use the
        # default build folder, otherwise the built object would be cached
        # and prevent recompilation
        _build_temp = self.build_temp
        self.build_temp = os.path.join(_build_temp, ext.name)
        _build_ext.build_extension(self, ext)
        self.build_temp = _build_temp

    def build_extensions(self):
        # check `cythonize` is available
        if isinstance(cythonize, ImportError):
            raise RuntimeError(
                "Cython is required to run `build_ext` command"
            ) from cythonize

        # remove universal compilation flags for OSX
        if platform.system() == "Darwin":
            _patch_osx_compiler(self.compiler)

        # use debug directives with Cython if building in debug mode
        cython_args = {
            "include_path": ["include"],
            "compiler_directives": {
                "cdivision": True,
                "nonecheck": False,
            },
            "compile_time_env": {
                "SYS_IMPLEMENTATION_NAME": sys.implementation.name,
                "SYS_VERSION_INFO_MAJOR": sys.version_info.major,
                "SYS_VERSION_INFO_MINOR": sys.version_info.minor,
                "SYS_VERSION_INFO_MICRO": sys.version_info.micro,
                "DEFAULT_BUFFER_SIZE": io.DEFAULT_BUFFER_SIZE,
                "TARGET_CPU": TARGET_CPU,
                "TARGET_SYSTEM": TARGET_SYSTEM,
                "AVX2_BUILD_SUPPORT": False,
                "NEON_BUILD_SUPPORT": False,
                "SSE2_BUILD_SUPPORT": False,
                "SSE4_BUILD_SUPPORT": False,
            },
        }
        if self.force:
            cython_args["force"] = True
        if self.debug:
            cython_args["annotate"] = True
            cython_args["compiler_directives"]["cdivision_warnings"] = True
            cython_args["compiler_directives"]["warn.undeclared"] = True
            cython_args["compiler_directives"]["warn.unreachable"] = True
            cython_args["compiler_directives"]["warn.maybe_uninitialized"] = True
            cython_args["compiler_directives"]["warn.unused"] = True
            cython_args["compiler_directives"]["warn.unused_arg"] = True
            cython_args["compiler_directives"]["warn.unused_result"] = True
            cython_args["compiler_directives"]["warn.multiple_declarators"] = True
        else:
            cython_args["compiler_directives"]["boundscheck"] = False
            cython_args["compiler_directives"]["wraparound"] = False

        # compile the C libraries
        if not self.distribution.have_run.get("build_clib", False):
            self._clib_cmd.run()

        # add the include dirs
        for ext in self.extensions:
            ext.include_dirs.append(self._clib_cmd.build_clib)

        # check if we can build platform-specific code
        if TARGET_CPU == "x86":
            if not self._simd_disabled["AVX2"] and self._check_avx2():
                cython_args["compile_time_env"]["AVX2_BUILD_SUPPORT"] = True
                self._simd_supported["AVX2"] = True
                self._simd_flags["AVX2"].extend(self._avx2_flags())
                self._simd_defines["AVX2"].append(("__AVX2__", 1))
            if not self._simd_disabled["SSE4"] and self._check_sse4():
                cython_args["compile_time_env"]["SSE4_BUILD_SUPPORT"] = True
                self._simd_supported["SSE4"] = True
                self._simd_flags["SSE4"].extend(self._sse4_flags())
                self._simd_defines["SSE4"].append(("__SSE4_1__", 1))
            if not self._simd_disabled["SSE2"] and self._check_sse2():
                cython_args["compile_time_env"]["SSE2_BUILD_SUPPORT"] = True
                self._simd_supported["SSE2"] = True
                self._simd_flags["SSE2"].extend(self._sse2_flags())
                self._simd_defines["SSE2"].append(("__SSE2__", 1))
        elif TARGET_CPU == "arm" or TARGET_CPU == "aarch64":
            if not self._simd_disabled["NEON"] and self._check_neon():
                cython_args["compile_time_env"]["NEON_BUILD_SUPPORT"] = True
                self._simd_supported["NEON"] = True
                self._simd_flags["NEON"].extend(self._neon_flags())
                self._simd_defines["NEON"].append(("__ARM_NEON__", 1))

        # filter out extensions missing required CPU features
        extensions = []
        for ext in self.extensions:
            if ext.requires is None:
                extensions.append(ext)
            elif self._simd_supported[ext.requires]:
                ext.extra_compile_args.extend(self._simd_flags[ext.requires])
                ext.extra_link_args.extend(self._simd_flags[ext.requires])
                extensions.append(ext)
        if not extensions:
            raise RuntimeError("Cannot build Opal for platform {}, no SIMD backend supported".format(MACHINE))

        # cythonize the extensions (retaining platform-specific sources)
        self.extensions = cythonize(extensions, **cython_args)

        # build the extensions as normal
        _build_ext.build_extensions(self)


class clean(_clean):
    """A `clean` that removes intermediate files created by Cython."""

    def run(self):

        source_dir = os.path.join(os.path.dirname(__file__), "pymuscle5")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c", "*.cpp"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                _eprint("removing {!r}".format(file))
                os.remove(file)

        _clean.run(self)


# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    ext_modules=[
        Extension(
            "pyopal._opal_neon",
            language="c++",
            requires="NEON",
            define_macros=[
                ("__ARM_NEON", 1),
            ],
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "_opal_neon.pyx"),
            ],
            include_dirs=[
                os.path.join("vendor", "opal", "src"),
                "pyopal",
                "include",
            ],
        ),
        Extension(
            "pyopal._opal_sse2",
            language="c++",
            requires="SSE2",
            define_macros=[
                ("__SSE2__", 1),
            ],
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "_opal_sse2.pyx"),
            ],
            include_dirs=[
                os.path.join("vendor", "opal", "src"),
                "pyopal",
                "include",
            ],
        ),
        Extension(
            "pyopal._opal_sse4",
            language="c++",
            requires="SSE4",
            define_macros=[
                ("__SSE4_1__", 1),
            ],
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "_opal_sse4.pyx"),
            ],
            include_dirs=[
                os.path.join("vendor", "opal", "src"),
                "pyopal",
                "include",
            ],
        ),
        Extension(
            "pyopal._opal_avx2",
            language="c++",
            requires="AVX2",
            define_macros=[
                ("__AVX2__", 1),
            ],
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "_opal_avx2.pyx"),
            ],
            include_dirs=[
                os.path.join("vendor", "opal", "src"),
                "pyopal",
                "include",
            ],
        ),
        Extension(
            "pyopal._opal",
            language="c++",
            sources=[
                os.path.join("vendor", "opal", "src", "ScoreMatrix.cpp"),
                os.path.join("pyopal", "_opal.pyx"),
            ],
            extra_compile_args=[],
            include_dirs=[
                os.path.join("vendor", "opal", "src"),
                "pyopal",
                "include",
            ],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "clean": clean,
    },
)
