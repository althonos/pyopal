# noqa: D104
from ._version import __version__

__author__ = "Martin Larralde <martin.larralde@embl.de>"
__license__ = "MIT"

from . import lib
from .lib import (
    Alphabet,
    Database,
    ScoreMatrix,
    ScoreResult,
    EndResult,
    FullResult,
)

__doc__ = lib.__doc__

# Small addition to the docstring: we want to show a link redirecting to the
# rendered version of the documentation, but this can only work when Python
# is running with docstrings enabled
if __doc__ is not None:
    __doc__ += """See Also:
    An online rendered version of the documentation for this version
    of the library on
    `Read The Docs <https://pyopal.readthedocs.io/en/v{}/>`_.

    """.format(
        __version__
    )
