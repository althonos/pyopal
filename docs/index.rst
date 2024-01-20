PyOpal |Stars|
==============

.. |Stars| image:: https://img.shields.io/github/stars/althonos/pyopal.svg?style=social&maxAge=3600&label=Star
   :target: https://github.com/althonos/pyopal/stargazers

`Cython <https://cython.org/>`_ *bindings and Python interface to* `Opal <https://github.com/Martinsos/opal>`_,
*a SIMD-accelerated pairwise aligner.*

|Actions| |Coverage| |PyPI| |Bioconda| |AUR| |Wheel| |Versions| |Implementations| |License| |Source| |Mirror| |Issues| |Docs| |Changelog| |Downloads|

.. |Actions| image:: https://img.shields.io/github/actions/workflow/status/althonos/pyopal/test.yml?branch=main&logo=github&style=flat-square&maxAge=300
   :target: https://github.com/althonos/pyopal/actions

.. |Coverage| image:: https://img.shields.io/codecov/c/gh/althonos/pyopal?style=flat-square&maxAge=600
   :target: https://codecov.io/gh/althonos/pyopal/

.. |PyPI| image:: https://img.shields.io/pypi/v/pyopal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.python.org/pypi/pyopal

.. |Bioconda| image:: https://img.shields.io/conda/vn/bioconda/pyopal?style=flat-square&maxAge=3600
   :target: https://anaconda.org/bioconda/pyopal

.. |AUR| image:: https://img.shields.io/aur/version/python-pyopal?logo=archlinux&style=flat-square&maxAge=3600
   :target: https://aur.archlinux.org/packages/python-pyopal

.. |Wheel| image:: https://img.shields.io/pypi/wheel/pyopal?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyopal/#files

.. |Versions| image:: https://img.shields.io/pypi/pyversions/pyopal.svg?style=flat-square&maxAge=3600
   :target: https://pypi.org/project/pyopal/#files

.. |Implementations| image:: https://img.shields.io/pypi/implementation/pyopal.svg?style=flat-square&maxAge=3600&label=impl
   :target: https://pypi.org/project/pyopal/#files

.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg?style=flat-square&maxAge=3600
   :target: https://choosealicense.com/licenses/mit/

.. |Source| image:: https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyopal/

.. |Mirror| image:: https://img.shields.io/badge/mirror-EMBL-009f4d?style=flat-square&maxAge=2678400
   :target: https://git.embl.de/larralde/pytrimapyopall/

.. |Issues| image:: https://img.shields.io/github/issues/althonos/pyopal.svg?style=flat-square&maxAge=600
   :target: https://github.com/althonos/pyopal/issues

.. |Docs| image:: https://img.shields.io/readthedocs/pyopal?style=flat-square&maxAge=3600
   :target: http://pyopal.readthedocs.io/en/stable/?badge=stable

.. |Changelog| image:: https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square
   :target: https://github.com/althonos/pyopal/blob/main/CHANGELOG.md

.. |Downloads| image:: https://img.shields.io/pypi/dm/pyopal?style=flat-square&color=303f9f&maxAge=86400&label=downloads
   :target: https://pepy.tech/project/pyopal


Overview
--------

PyOpal is a Python module that provides bindings to `Opal <https://github.com/Martinsos/opal>`_ 
using `Cython <https://cython.org/>`_. It directly interacts with the trimAl
internals, which has the following advantages:

- **no binary dependency**: PyOpal is distributed as a Python package, so you
  can add it as a dependency to your project, and stop worrying about the
  Opal binary being present on the end-user machine.
- **no intermediate files**: Everything happens in memory, in a Python object
  you control, so you don't have to invoke the Opal CLI using a sub-process
  and temporary files.
- **better portability**: Opal uses SIMD to accelerate alignment scoring, but
  doesn't support dynamic dispatch, so it has to be compiled on the local
  machine to be able to use the full capabilities of the local CPU. PyOpal
  ships several versions of Opal instead, each compiled with different target
  features, and selects the best one for the local platform at runtime.
- **wider platform support**: The Opal code has been backported to work on SSE2
  rather than SSE4.1, allowing PyOpal to run on older x86 CPUs (all x86 CPUs
  support it since 2003). In addition, Armv7 and Aarch64 CPUs are also 
  supported if they implement NEON extensions. Finally, the C++ code of Opal
  has been modified to compile on Windows.


Setup
-----

Run ``pip install pyopal`` in a shell to download the latest release from PyPi, 
or have a look at the :doc:`Installation page <install>` to find other ways 
to install ``pyopal``.


Library
-------

.. toctree::
   :maxdepth: 2

   Installation <install>
   Examples <examples/index>
   Contributing <contributing>
   API Reference <api/index>
   Changelog <changes>


License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
Opal was developed by `Martin Šošić <https://github.com/Martinsos>`_ and is distributed under the
terms of the MIT License as well. 

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `Opal`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
