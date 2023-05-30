# noqa: D104

from . import (
    test_database,
    test_doctest,
    test_search,
)


def load_tests(loader, suite, pattern):
    test_doctest.load_tests(loader, suite, pattern)
    suite.addTests(loader.loadTestsFromModule(test_database))
    suite.addTests(loader.loadTestsFromModule(test_search))
    return suite
