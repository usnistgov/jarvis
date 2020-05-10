"""
Helper functions for JSON files
"""

import json


def loadjson(filename=""):
    """
    Helper function to load a json file
    Serialization would be available soon
    """
    f = open(filename, "r")
    d = json.load(f)
    f.close()
    return d


def dumpjson(data=[], filename=""):
    """
    Helper function to write a json file
    Serialization would be available soon
    """
    f = open(filename, "w")
    f.write(json.dumps(data))
    f.close()
