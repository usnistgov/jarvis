import json


def loadjson(filename=""):
    f = open(filename, "r")
    d = json.load(f)
    f.close()
    return d


def dumpjson(data=[], filename="", indent=0):
    f = open(filename, "w")
    f.write(json.dumps(data, indent=indent))
    f.close()
