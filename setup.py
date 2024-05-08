import configparser
import functools
import glob
import itertools
import io
import json
import multiprocessing.pool
import os
import platform
import re
import setuptools
import setuptools.extension
import subprocess
import string
import sys
import sysconfig
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
        self.simd = kwargs.pop("simd", None)
        super().__init__(*args, **kwargs)

class ExtensionTemplate(Extension):
    def __init__(self, *args, **kwargs):
        self.templates = kwargs.pop("templates", {})
        super().__init__(*args, **kwargs)


# --- Utils ------------------------------------------------------------------

_HEADER_PATTERN = re.compile("^@@ -(\d+),?(\d+)? \+(\d+),?(\d+)? @@$")


def _eprint(*args, **kwargs):
    print(*args, **kwargs, file=sys.stderr)


def _patch_osx_compiler(compiler, machine):
    # On newer OSX, Python has been compiled as a universal binary, so
    # it will attempt to pass universal binary flags when building the
    # extension. This will not work because the code makes use of SSE2.
    for tool in ("compiler", "compiler_so", "linker_so"):
        flags = getattr(compiler, tool)
        i = next(
            (
                i
                for i in range(1, len(flags))
                if flags[i - 1] == "-arch" and flags[i] != machine
            ),
            None,
        )
        if i is not None:
            flags.pop(i)
            flags.pop(i - 1)


def _detect_target_machine(platform):
    if platform == "win32":
        return "x86"
    return platform.rsplit("-", 1)[-1]


def _detect_target_cpu(platform):
    machine = _detect_target_machine(platform)
    if re.match("^mips", machine):
        return "mips"
    elif re.match("^(aarch64|arm64)$", machine):
        return "aarch64"
    elif re.match("^arm", machine):
        return "arm"
    elif re.match("(x86_64)|AMD64|amd64", machine):
        return "x86_64"
    elif re.match("(x86)|(^i.86$)", machine):
        return "x86"
    elif re.match("^(powerpc|ppc)", machine):
        return "ppc"
    return None


def _detect_target_system(platform):
    if platform.startswith("win"):
        return "windows"
    elif platform.startswith("macos"):
        return "macos"
    elif platform.startswith("linux"):
        return "linux_or_android"
    elif platform.startswith("freebsd"):
        return "freebsd"
    return None


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
        self.target_machine = None
        self.target_system = None
        self.target_cpu = None

    def finalize_options(self):
        _build_ext.finalize_options(self)
        # check platform
        if self.plat_name is None:
            self.plat_name = sysconfig.get_platform()
        # detect platform options
        self.target_machine = _detect_target_machine(self.plat_name)
        self.target_system = _detect_target_system(self.plat_name)
        self.target_cpu = _detect_target_cpu(self.plat_name)
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
        except (subprocess.CalledProcessError, OSError):
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
        return ["-mfpu=neon"] if self.target_cpu == "arm" else []

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

    # --- Template code ---

    def substitute_template(self, template_file, output_file, metadata):
        with open(template_file) as f:
            template = string.Template(f.read())
        with open(output_file, "w") as dst:
            dst.write(template.substitute(metadata))

    def generate_extension(self, ext):
        metadata = dict(SIMD=ext.simd)
        for out_file, template_file in ext.templates.items():
            self.make_file([template_file], out_file, self.substitute_template, (template_file, out_file, metadata))

    # --- Build code ---

    def build_extension(self, ext):
        # show the compiler being used
        _eprint("building", ext.name, "with", self.compiler.compiler_type, "compiler for platform", self.plat_name)

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
                "-faligned-new",
                "-std=c++17",
                "-Wno-unused-variable",
                "-Wno-maybe-uninitialized",
                "-Wno-return-type"
            ])
        elif self.compiler.compiler_type == "msvc":
            ext.extra_compile_args.append("/std:c++17")

        # add Windows flags
        if self.compiler.compiler_type == "msvc":
            ext.define_macros.append(("WIN32", 1))

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
            _patch_osx_compiler(self.compiler, self.target_machine)

        # generate files from templates:
        for i, ext in enumerate(self.extensions):
            if isinstance(ext, ExtensionTemplate):
                self.generate_extension(ext)

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
                "TARGET_CPU": self.target_cpu,
                "TARGET_SYSTEM": self.target_system,
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

        # check if `PyInterpreterState_GetID` is defined
        if self._check_getid():
            self.compiler.define_macro("HAS_PYINTERPRETERSTATE_GETID", 1)

        # check if we can build platform-specific code
        if self.target_cpu == "x86" or self.target_cpu == "x86_64":
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
        elif self.target_cpu == "arm" or self.target_cpu == "aarch64":
            if not self._simd_disabled["NEON"] and self._check_neon():
                cython_args["compile_time_env"]["NEON_BUILD_SUPPORT"] = True
                self._simd_supported["NEON"] = True
                self._simd_flags["NEON"].extend(self._neon_flags())
                self._simd_defines["NEON"].append(("__ARM_NEON__", 1))

        # filter out extensions missing required CPU features
        extensions = []
        for ext in self.extensions:
            if ext.simd is None:
                extensions.append(ext)
            elif self._simd_supported[ext.simd]:
                ext.extra_compile_args.extend(self._simd_flags[ext.simd])
                ext.extra_link_args.extend(self._simd_flags[ext.simd])
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

        source_dir = os.path.join(os.path.dirname(__file__), "pyopal")

        patterns = ["*.html"]
        if self.all:
            patterns.extend(["*.so", "*.c", "*.cpp"])

        for pattern in patterns:
            for file in glob.glob(os.path.join(source_dir, pattern)):
                _eprint("removing {!r}".format(file))
                os.remove(file)

        for ext in self.distribution.ext_modules:
            if isinstance(ext, ExtensionTemplate):
                for source_file in ext.templates:
                    if os.path.exists(source_file):
                        _eprint("removing {!r}".format(source_file))
                        os.remove(source_file)

            for source_file in ext.sources:
                if source_file.endswith(".pyx"):
                    ext = ".cpp" if ext.language == "c++" else ".c"
                    c_file = source_file.replace(".pyx", ext)
                    if os.path.exists(c_file):
                        _eprint("removing {!r}".format(c_file))
                        os.remove(c_file)

        _clean.run(self)


# --- Setup ---------------------------------------------------------------------

setuptools.setup(
    ext_modules=[
        ExtensionTemplate(
            "pyopal.platform.neon",
            language="c++",
            simd="NEON",
            define_macros=[("__ARM_NEON", 1)],
            include_dirs=["pyopal"],
            templates={
                os.path.join("pyopal", "platform", "neon.pxd"): os.path.join("pyopal", "platform", "pxd.in"),
                os.path.join("pyopal", "platform", "neon.pyx"): os.path.join("pyopal", "platform", "pyx.in"),
            },
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "platform", "neon.pyx"),
            ],
        ),
        ExtensionTemplate(
            "pyopal.platform.sse2",
            language="c++",
            simd="SSE2",
            define_macros=[("__SSE2__", 1)],
            include_dirs=["pyopal"],
            templates={
                os.path.join("pyopal", "platform", "sse2.pxd"): os.path.join("pyopal", "platform", "pxd.in"),
                os.path.join("pyopal", "platform", "sse2.pyx"): os.path.join("pyopal", "platform", "pyx.in"),
            },
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "platform", "sse2.pyx"),
            ],
        ),
        ExtensionTemplate(
            "pyopal.platform.sse4",
            language="c++",
            simd="SSE4",
            define_macros=[("__SSE4_1__", 1)],
            include_dirs=["pyopal"],
            templates={
                os.path.join("pyopal", "platform", "sse4.pxd"): os.path.join("pyopal", "platform", "pxd.in"),
                os.path.join("pyopal", "platform", "sse4.pyx"): os.path.join("pyopal", "platform", "pyx.in"),
            },
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "platform", "sse4.pyx"),
            ],
        ),
        ExtensionTemplate(
            "pyopal.platform.avx2",
            language="c++",
            simd="AVX2",
            include_dirs=["pyopal"],
            define_macros=[("__AVX2__", 1)],
            templates={
                os.path.join("pyopal", "platform", "avx2.pxd"): os.path.join("pyopal", "platform", "pxd.in"),
                os.path.join("pyopal", "platform", "avx2.pyx"): os.path.join("pyopal", "platform", "pyx.in"),
            },
            sources=[
                os.path.join("vendor", "opal", "src", "opal.cpp"),
                os.path.join("pyopal", "platform", "avx2.pyx"),
            ],
        ),
        Extension(
            "pyopal.lib",
            language="c++",
            sources=[os.path.join("pyopal", "lib.pyx")],
            include_dirs=["pyopal"],
        ),
    ],
    cmdclass={
        "sdist": sdist,
        "build_ext": build_ext,
        "clean": clean,
    },
)
