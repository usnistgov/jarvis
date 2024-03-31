from jarvis.db.restapi import jarvisdft_optimade


def test_optimade():
    x = jarvisdft_optimade()
    print(x)
    print(len(x))
    assert len(x) > 1
