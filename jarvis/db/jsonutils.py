"""Helper functions for JSON files."""

import json


def loadjson(filename=""):
    """Provide helper function to load a json file."""
    f = open(filename, "r")
    d = json.load(f)
    f.close()
    return d


def dumpjson(data=[], filename=""):
    """Provide helper function to write a json file."""
    f = open(filename, "w")
    f.write(json.dumps(data))
    f.close()
