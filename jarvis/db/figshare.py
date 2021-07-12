"""
Downloads files from Figshare.

Main page: https://figshare.com/authors/Kamal_Choudhary/4445539
"""

import zipfile
import tempfile
import os
import numpy as np
import io
import json
import requests
from jarvis.db.jsonutils import loadjson
from jarvis.io.vasp.outputs import Vasprun
from jarvis.io.vasp.inputs import Poscar
from jarvis.io.wannier.outputs import WannierHam
from tqdm import tqdm
import matplotlib.image as mpimg
from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
import matplotlib.pyplot as plt

plt.switch_backend("agg")


def get_db_info():
    """Get DB info."""
    db_info = {
        "dft_2d": [
            "https://ndownloader.figshare.com/files/26808917",
            "d2-3-12-2021.json",
            "Obtaining 2D dataset 1.1k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.6815705
        "dft_3d": [
            "https://ndownloader.figshare.com/files/28034994",
            "d3-5-16-2021.json",
            "Obtaining 3D dataset 48k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.6815699
        "cfid_3d": [
            "https://ndownloader.figshare.com/files/26808914",
            "d3-3-12-2021.json",
            "Obtaining 3D dataset 43k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        "jff": [
            "https://ndownloader.figshare.com/files/26809760",
            "jff-3-12-2021.json",
            "Obtaining JARVIS-FF 2k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.14213522
        "mp_3d_2020": [
            "https://ndownloader.figshare.com/files/26791259",
            "all_mp.json",
            "Obtaining Materials Project-3D CFID dataset 127k...",
            "https://doi.org/10.1063/1.4812323",
        ],
        "megnet": [
            "https://ndownloader.figshare.com/files/26724977",
            "megnet.json",
            "Obtaining MEGNET-3D CFID dataset 69k...",
            "https://pubs.acs.org/doi/10.1021/acs.chemmater.9b01294",
        ],
        # https://doi.org/10.6084/m9.figshare.14177630
        "megnet2": [
            "https://ndownloader.figshare.com/files/28332741",
            "megnet-mp-2019-04-01.json",
            "Obtaining MEGNET-3D CFID dataset 133k...",
            "https://pubs.acs.org/doi/10.1021/acs.chemmater.9b01294",
        ],
        # https://doi.org/10.6084/m9.figshare.14745435
        "edos_pdos": [
            "https://ndownloader.figshare.com/files/28501764",
            "edos-up_pdos-elast_interp-6-19-2021.json",
            "Interpolated electronic total dos spin-up dataset 48k...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.14745327
        "mp_3d": [
            "https://ndownloader.figshare.com/files/24979850",
            "CFID_mp_desc_data_84k.json",
            "Obtaining Materials Project-3D CFID dataset 84k...",
            "https://doi.org/10.1063/1.4812323",
        ],
        # https://doi.org/10.6084/m9.figshare.13054247
        "oqmd_3d": [
            "https://ndownloader.figshare.com/files/24981170",
            "CFID_OQMD_460k.json",
            "Obtaining OQMD-3D CFID dataset 460k...",
            "https://www.nature.com/articles/npjcompumats201510",
        ],
        # https://doi.org/10.6084/m9.figshare.13055333
        "oqmd_3d_no_cfid": [
            "https://ndownloader.figshare.com/files/26790182",
            "all_oqmd.json",
            "Obtaining OQMD-3D  dataset 800k...",
            "https://www.nature.com/articles/npjcompumats201510",
        ],
        # https://doi.org/10.6084/m9.figshare.14206169
        "twod_matpd": [
            "https://ndownloader.figshare.com/files/26789006",
            "twodmatpd.json",
            "Obtaining 2DMatPedia dataset 6k...",
            "https://www.nature.com/articles/s41597-019-0097-3",
        ],
        # https://doi.org/10.6084/m9.figshare.14205083
        "polymer_genome": [
            "https://ndownloader.figshare.com/files/26809907",
            "pgnome.json",
            "Obtaining Polymer genome 1k...",
            "https://www.nature.com/articles/sdata201612",
        ],
        # https://doi.org/10.6084/m9.figshare.14213603
        "qm9_std_jctc": [
            "https://ndownloader.figshare.com/files/28715319",
            "qm9_std_jctc.json",
            "Obtaining QM9 standardized dataset 130k,"
            + "From https://doi.org/10.1021/acs.jctc.7b00577,+",
            "https://www.nature.com/articles/sdata201422",
        ],
        "qm9_dgl": [
            "https://ndownloader.figshare.com/files/28541196",
            "qm9_dgl.json",
            "Obtaining QM9 dataset 130k, from DGL...",
            "https://www.nature.com/articles/sdata201422",
        ],
        # https://doi.org/10.6084/m9.figshare.14827584
        # Use qm9_std_jctc instaed
        "cod": [
            "https://ndownloader.figshare.com/files/28715301",
            "cod_db.json",
            "Obtaining COD dataset 431k",
            "https://doi.org/10.1107/S1600576720016532",
        ],
        # https://doi.org/10.6084/m9.figshare.14912820.v1
        "qm9": [
            "https://ndownloader.figshare.com/files/27627596",
            "qm9_data_cfid.json",
            "Obtaining QM9 dataset 134k...",
            "https://www.nature.com/articles/sdata201422",
        ],
        # Use qm9_std_jctc instaed
        "qe_tb": [
            "https://ndownloader.figshare.com/files/TODO",
            "jqe_tb_folder.json",
            "Obtaining QETB dataset 860k...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        "omdb": [
            "https://ndownloader.figshare.com/files/28501761",
            "omdbv1.json",
            "Obtaining OMDB dataset 12.5k...",
            "https://doi.org/10.1002/qute.201900023",
        ],
        # https://doi.org/10.6084/m9.figshare.14812050
        "qmof": [
            "https://ndownloader.figshare.com/files/28501740",
            "qmof_db.json",
            "Obtaining QMOF dataset 18k...",
            "https://www.cell.com/matter/fulltext/S2590-2385(21)00070-9",
        ],
        # https://doi.org/10.6084/m9.figshare.14812044
        "c2db": [
            "https://ndownloader.figshare.com/files/28682010",
            "c2db_atoms.json",
            "Obtaining C2DB dataset 3.5k...",
            "https://iopscience.iop.org/article/10.1088/2053-1583/aacfc1",
        ],
        # https://figshare.com/account/projects/100325/articles/14960157
        "hpov": [
            "https://ndownloader.figshare.com/files/28814184",
            "hopv_15.json",
            "Obtaining HPOV15 dataset 4.5k...",
            "https://www.nature.com/articles/sdata201686",
        ],
        # https://figshare.com/account/projects/100325/articles/14962356
        "pdbbind": [
            "https://ndownloader.figshare.com/files/28816368",
            "pdbbind_2015.json",
            "Obtaining PDBBind dataset 11k...",
            "https://doi.org/10.1093/bioinformatics/btu626",
        ],
        # https://doi.org/10.6084/m9.figshare.14812038
        "aflow2": [
            "https://ndownloader.figshare.com/files/25453265",
            "CFID_AFLOW2.json",
            "Obtaining AFLOW-2 CFID dataset 400k...",
            "https://doi.org/10.1016/j.commatsci.2012.02.005",
        ],
        # https://doi.org/10.6084/m9.figshare.13215308
        "arXiv": [
            "https://ndownloader.figshare.com/files/26804795",
            "arXivdataset.json",
            "Obtaining arXiv dataset 1.8 million...",
            "https://www.kaggle.com/Cornell-University/arxiv",
        ],
        # https://doi.org/10.6084/m9.figshare.14211860
        "cord19": [
            "https://ndownloader.figshare.com/files/26804798",
            "cord19.json",
            "Obtaining CORD19 dataset 223k...",
            "https://github.com/usnistgov/cord19-cdcs-nist",
        ],
        # https://doi.org/10.6084/m9.figshare.14211857
        "raw_files": [
            "https://ndownloader.figshare.com/files/25295732",
            "figshare_data-10-28-2020.json",
            "Obtaining raw io files 145k...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
    }
    return db_info


# Format: download_link, filename, message, reference
# Figshare link: https://figshare.com/account/home#/projects/100325


def get_stm_2d_dataset():
    """Get 2D STM image dataset."""
    # Ref: https://www.nature.com/articles/s41597-021-00824-y
    link_1 = "https://ndownloader.figshare.com/files/21884952"
    r_jpg = requests.get(link_1)
    z = zipfile.ZipFile(io.BytesIO(r_jpg.content))
    link_2 = "https://ndownloader.figshare.com/files/21893379"
    r_json = requests.get(link_2).content
    latts = json.loads(r_json)
    namelist = z.namelist()
    pos_bias = []
    neg_bias = []

    print("Obtaining 2D STM dataset ...")
    for i in namelist:
        img_str = z.read(i)
        values = mpimg.imread(io.BytesIO((img_str)), format="jpg")
        # img=Image(values=values)
        jid = i.split("/")[-1].split("_")[0]
        bias = i.split("/")[-1].split("_")[1].split(".jpg")[0]
        lat_system = latts[jid]
        if bias == "pos":
            info = {}
            info["jid"] = jid
            info["image_values"] = values
            info["lat_type"] = lat_system
            pos_bias.append(info)
        if bias == "neg":
            info = {}
            info["jid"] = jid
            info["image_values"] = values
            info["lat_type"] = lat_system
            neg_bias.append(info)
    return pos_bias, neg_bias


def get_request_data(
    js_tag="jdft_2d-4-26-2020.json",
    url="https://ndownloader.figshare.com/files/22471019",
):
    """Get data with progress bar."""
    zfile = js_tag + ".zip"
    path = str(os.path.join(os.path.dirname(__file__), zfile))
    # path = str(os.path.join(os.path.dirname(__file__), js_tag))
    if not os.path.isfile(path):
        # zfile = str(os.path.join(os.path.dirname(__file__), "tmp.zip"))
        response = requests.get(url, stream=True)
        total_size_in_bytes = int(response.headers.get("content-length", 0))
        block_size = 1024  # 1 Kibibyte
        progress_bar = tqdm(
            total=total_size_in_bytes, unit="iB", unit_scale=True
        )
        with open(path, "wb") as file:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))
                file.write(data)
        progress_bar.close()
        # f = open(zfile, "wb")
        # f.write(r.content)
        # f.close()
    print("Loading the zipfile...")
    data = json.loads(zipfile.ZipFile(path).read(js_tag))
    print("Loading completed.")

    #    with zipfile.ZipFile(zfile, "r") as zipObj:
    #        # zipObj.extract(path)
    #        zipObj.extractall(os.path.join(os.path.dirname(__file__)))
    #    os.remove(zfile)
    # data = loadjson(path)
    return data


def data(dataset="dft_2d"):
    """Provide main function to download datasets."""
    db_info = get_db_info()
    if dataset not in list(db_info.keys()):
        raise ValueError("Check DB name options.")
    url = db_info[dataset][0]
    js_tag = db_info[dataset][1]
    message = db_info[dataset][2]
    print(message)
    message = "Reference:" + db_info[dataset][3]
    print(message)
    # r = requests.get(url)
    # z = zipfile.ZipFile(io.BytesIO(r.content))
    # data = json.loads(z.read(js_tag).decode("utf-8"))

    # r = requests.get(url)
    # z = zipfile.ZipFile(io.BytesIO(r.content))
    # wdat = z.read(js_tag).decode("utf-8")
    # fd, path = tempfile.mkstemp()
    # with os.fdopen(fd, "w") as tmp:
    #    tmp.write(wdat)
    # data = loadjson(path)

    # path = str(os.path.join(os.path.dirname(__file__), js_tag))
    # if not os.path.isfile(path):
    #    zfile = str(os.path.join(os.path.dirname(__file__), "tmp.zip"))
    #    r = requests.get(url)
    #    f = open(zfile, "wb")
    #    f.write(r.content)
    #    f.close()

    #    with zipfile.ZipFile(zfile, "r") as zipObj:
    #        # zipObj.extract(path)
    #        zipObj.extractall(os.path.join(os.path.dirname(__file__)))
    #    os.remove(zfile)
    # data = loadjson(path)
    dat = get_request_data(js_tag=js_tag, url=url)
    return dat


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


def make_stm_from_prev_parchg(
    jid="JVASP-667", bias="Negative", filename="stm_image.png", min_size=10
):
    """Make STM images from previously calculated PARVHG files for 2D."""
    fls = data("raw_files")
    for i in fls["STM"]:
        zip_name = jid + "_" + bias + ".zip"
        if i["name"] == zip_name:
            zip_file_url = i["download_url"]
            r = requests.get(zip_file_url)
            z = zipfile.ZipFile(io.BytesIO(r.content))
            pchg = z.read("PARCHG").decode("utf-8")
            fd, path = tempfile.mkstemp()
            with os.fdopen(fd, "w") as tmp:
                tmp.write(pchg)
            TH_STM = TersoffHamannSTM(
                chg_name=path, min_size=min_size, zcut=None
            )
            t_height = TH_STM.constant_height(filename=filename)
            print("t_height", t_height)
            return i


def get_wann_electron(jid="JVASP-816"):
    """Download electron WTBH if available."""
    w = ""
    ef = ""
    fls = data("raw_files")
    for i in fls["WANN"]:
        if i["name"].split(".zip")[0] == jid:
            r = requests.get(i["download_url"])
            z = zipfile.ZipFile(io.BytesIO(r.content))
            wdat = z.read("wannier90_hr.dat").decode("utf-8")
            js_file = jid + ".json"
            js = z.read(js_file).decode("utf-8")
            fd, path = tempfile.mkstemp()
            with os.fdopen(fd, "w") as tmp:
                tmp.write(wdat)
            w = WannierHam(path)
            fd, path = tempfile.mkstemp()
            with os.fdopen(fd, "w") as tmp:
                tmp.write(js)
            d = loadjson(path)
            ef = d["info_mesh"]["efermi"]
            fd, path = tempfile.mkstemp()
            pos = z.read("POSCAR").decode("utf-8")
            with os.fdopen(fd, "w") as tmp:
                tmp.write(pos)
            atoms = Poscar.from_file(path).atoms
    return w, ef, atoms


def get_wann_phonon(jid="JVASP-1002", factor=15.633302):
    """Download phonon WTBH if available."""
    # Requires phonopy
    from jarvis.io.phonopy.outputs import get_phonon_tb

    fls = data("raw_files")

    for i in fls["FD-ELAST"]:
        if isinstance(i, dict):
            if i["name"].split(".zip")[0] == jid:
                r = requests.get(i["download_url"])
                z = zipfile.ZipFile(io.BytesIO(r.content))
                vrun_path = z.read("vasprun.xml").decode("utf-8")
                fd, path = tempfile.mkstemp()
                with os.fdopen(fd, "w") as tmp:
                    tmp.write(vrun_path)
                vrun = Vasprun(path)
                fc = vrun.phonon_data()["force_constants"]
                atoms = vrun.all_structures[0]
                # print(atoms)
                # atoms = Atoms.from_poscar(pos)
                # print(atoms)
                fd, path = tempfile.mkstemp()
                get_phonon_tb(fc=fc, atoms=atoms, out_file=path, factor=factor)
                # cvn = Spacegroup3D(atoms).conventional_standard_structure
                w = WannierHam(path)
                return w, atoms


def get_hk_tb(k=np.array([0, 0, 0]), w=[]):
    """Get Wannier TB Hamiltonian."""
    nr = w.R.shape[0]
    hk = np.zeros((w.nwan, w.nwan), dtype=complex)
    kmat = np.tile(k, (nr, 1))
    exp_ikr = np.exp(1.0j * 2 * np.pi * np.sum(kmat * w.R, 1))
    temp = np.zeros(w.nwan ** 2, dtype=complex)
    for i in range(nr):
        temp += exp_ikr[i] * w.HR[i, :]
    hk = np.reshape(temp, (w.nwan, w.nwan))
    hk = (hk + hk.T.conj()) / 2.0
    return hk


"""
QM9 xyz file
>>> def get_val(filename=''):
        f=open(filename,'r')
        lines=f.read().splitlines()
        f.close()
        info={}
        for i in range(len(attr_index)):
           info[attr_index[i]]=float(lines[1].split()[i+2])
        return info
"""

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
