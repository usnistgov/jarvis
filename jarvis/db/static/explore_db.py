from monty.json import MontyEncoder, MontyDecoder
from monty.serialization import loadfn, dumpfn
import os, requests

def get_ff_eneleast():
    jff1 = str(os.path.join(os.path.dirname(__file__), "jff1.json"))
    if not os.path.isfile(jff1):
        r = requests.get("https://ndownloader.figshare.com/files/10307139")
        f = open(jff1, "wb")
        f.write(r.content)
        f.close()
    data_ff1 = loadfn(jff1, cls=MontyDecoder)
    return data_ff1

def get_3d_dataset():
    j3d = str(os.path.join(os.path.dirname(__file__), "j3d.json"))
    if not os.path.isfile(j3d):
        r = requests.get("https://ndownloader.figshare.com/files/15475826")
        f = open(j3d, "wb")
        f.write(r.content)
        f.close()
    data_3d = loadfn(j3d, cls=MontyDecoder)
    return data_3d


def get_2d_dataset():
    j2d = str(os.path.join(os.path.dirname(__file__), "j2d.json"))
    if not os.path.isfile(j2d):
        r = requests.get("https://ndownloader.figshare.com/files/15475985")
        f = open(j2d, "wb")
        f.write(r.content)
    data_2d = loadfn(j2d, cls=MontyDecoder)
    return data_2d


def get_ml_dataset():
    jml = str(os.path.join(os.path.dirname(__file__), "jml.json"))
    if not os.path.isfile(jml):
        r = requests.get("https://ndownloader.figshare.com/files/12533579")
        f = open(jml, "wb")
        f.write(r.content)
    data_ml = loadfn(jml, cls=MontyDecoder)
    return data_ml


if __name__ == "__main__":

    data_2d = get_2d_dataset()
    print(len(data_2d))
    data_3d = get_3d_dataset()
    print(len(data_3d))
    data_ml = get_ml_dataset()
    print(len(data_ml))
    data_ff = get_ff_eneleast()
    print (len(data_ff))
