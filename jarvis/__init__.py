"""Version number."""
__version__ = "2021.12.16"

import os


def test(*args):

    import pytest

    path = os.path.join(os.path.split(__file__)[0], "tests")
    pytest.main(args=[path] + list(args))
