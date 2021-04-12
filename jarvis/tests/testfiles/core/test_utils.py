from jarvis.core.utils import (
    update_dict,
    recast_array_on_uniq_array_elements,
    stringdict_to_xml,
    array_to_string,
    check_url_exists,
)


def test_utils():
    info = {"a": "b", "c": "d"}
    sd = stringdict_to_xml(info)
    sd = stringdict_to_xml(info, enforce_string=True)
    ar = [1, 2, 3, 4, 5]
    sarr = array_to_string(ar)
    x = {"x": 1, "y": 2, "z": 3}
    y = {"m": 1, "n": 2, "o": 3}
    z = update_dict(x, y)
