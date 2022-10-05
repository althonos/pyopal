# noqa: D104

from . import (
    test_doctest,
    test_search,
)

def load_tests(loader, suite, pattern):
    suite.addTests(loader.loadTestsFromModule(test_doctest))
    suite.addTests(loader.loadTestsFromModule(test_search))
    return suite
