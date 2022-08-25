"""Version number."""
__version__ = "2022.07.17"

import os


def test(*args):
    """Run pytest in the base of jarvis."""
    import pytest

    path = os.path.join(os.path.split(__file__)[0], "tests")
    pytest.main(args=[path] + list(args))
