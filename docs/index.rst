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
   :target: https://git.embl.de/larralde/pyopal/

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
using `Cython <https://cython.org/>`_. It directly interacts with the Opal
internals, which has the following advantages:

.. grid:: 1 2 3 3
   :gutter: 1

   .. grid-item-card:: :fas:`battery-full` Batteries-included

      Just add ``pyopal`` as a ``pip`` or ``conda`` dependency, no need
      for the Opal binary or any external dependency.

   .. grid-item-card:: :fas:`screwdriver-wrench` Flexible

      Directly pass sequences to process as Python `str` objects, no 
      need for intermediate files.

   .. grid-item-card:: :fas:`dolly` Portable

      Use the best available SIMD implementation for your current platform
      without having to recompile Opal.

   .. grid-item-card:: :fas:`computer` Wider support

      The Opal code has been backported to work on SSE2 rather than SSE4.1, 
      allowing PyOpal to run on older x86 CPUs (all x86 CPUs support it 
      since 2003).

   .. grid-item-card:: :fas:`microchip` Arm-compatible

      Aarch64 CPUs (and Armv7 when they implement NEON extensions) are 
      supported in addition to x86-64 in the original.

   .. grid-item-card:: :fas:`arrow-pointer` Windows-compatible

      The C++ code of Opal has been modified to compile on Windows.


Setup
-----

Run ``pip install pyopal`` in a shell to download the latest release from PyPi, 
or have a look at the :doc:`Installation page <guide/install>` to find other ways 
to install ``pyopal``.


Library
-------

.. toctree::
   :maxdepth: 2

   User Guide <guide/index>
   Examples <examples/index>
   API Reference <api/index>


Related Projects
----------------

The following Python libraries may be of interest for bioinformaticians.

.. include:: related.rst


License
-------

This library is provided under the `MIT License <https://choosealicense.com/licenses/mit/>`_.
Opal was developed by `Martin Šošić <https://github.com/Martinsos>`_ and is distributed under the
terms of the MIT License as well. See the :doc:`Copyright Notice <guide/copyright>` section
for the full license.

*This project is in no way not affiliated, sponsored, or otherwise endorsed by
the original* `Opal`_ *authors. It was developed by* `Martin Larralde <https://github.com/althonos>`_ *during his
PhD project at the* `European Molecular Biology Laboratory <https://www.embl.de/>`_
*in the* `Zeller team <https://github.com/zellerlab>`_.
