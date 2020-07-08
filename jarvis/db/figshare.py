"""
Downloads files from Figshare.

Main page: https://figshare.com/authors/Kamal_Choudhary/4445539
"""

import zipfile
import os
import requests
from jarvis.db.jsonutils import loadjson


def datasets(dataset=""):
    """Get collection of dataset names and URLs."""
    if dataset == "dft_2d":
        url = "https://ndownloader.figshare.com/files/22471019"
        js_tag = "jdft_2d-4-26-2020.json"
        print("Obtaining 2D dataset ...")
    elif dataset == "dft_3d":
        url = "https://ndownloader.figshare.com/files/22471022"
        js_tag = "jdft_3d-4-26-2020.json"
        print("Obtaining 3D dataset ...")

    elif dataset == "cfid_3d":
        url = "https://ndownloader.figshare.com/files/22470818"
        js_tag = "jml_3d-4-26-2020.json"
        print("Obtaining 3D CFID dataset ...")
    else:
        ValueError("Dataset doesnt exist", dataset)
    return url, js_tag


def data(dataset="dft_2d"):
    """Provide main function to download datasets."""
    url, js_tag = datasets(dataset)
    path = str(os.path.join(os.path.dirname(__file__), js_tag))
    if not os.path.isfile(path):
        zfile = str(os.path.join(os.path.dirname(__file__), "tmp.zip"))
        r = requests.get(url)
        f = open(zfile, "wb")
        f.write(r.content)
        f.close()

        with zipfile.ZipFile(zfile, "r") as zipObj:
            # zipObj.extract(path)
            zipObj.extractall(os.path.join(os.path.dirname(__file__)))
        os.remove(zfile)
    data = loadjson(path)
    return data


def get_jid_data(jid="JVASP-667", dataset="dft_2d"):
    """Get info for a jid and dataset."""
    d = data(dataset)
    for i in d:
        if i["jid"] == jid:
            return i


def get_ff_eneleast():
    """Get JARVIS-FF related data."""
    jff1 = str(os.path.join(os.path.dirname(__file__), "jff1.json"))
    if not os.path.isfile(jff1):
        r = requests.get("https://ndownloader.figshare.com/files/10307139")
        f = open(jff1, "wb")
        f.write(r.content)
        f.close()
    data_ff1 = loadjson(jff1)
    return data_ff1


"""
if __name__ == "__main__":

    data_2d = data(dataset='dft_2d')
    print('2d',len(data_2d))
    data_3d = data(dataset='dft_3d')
    print('3d',len(data_3d))
    data_ml = data(dataset='cfid_3d')
    print('cfid3d',len(data_ml))
    data_ff = get_ff_eneleast()
    print ('ff',len(data_ff))
"""
