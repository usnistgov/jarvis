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
from tqdm import tqdm
import matplotlib.image as mpimg

# from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM
# from jarvis.io.wannier.outputs import WannierHam
# from jarvis.io.vasp.outputs import Vasprun
# from jarvis.io.vasp.inputs import Poscar

# import matplotlib.pyplot as plt
# plt.switch_backend("agg")


def get_db_info():
    """Get DB info."""
    db_info = {
        # https://doi.org/10.6084/m9.figshare.6815705
        "dft_2d": [
            "https://ndownloader.figshare.com/files/38521268",
            "d2-12-12-2022.json",
            "Obtaining 2D dataset 1.1k ...",
            "https://www.nature.com/articles/s41524-020-00440-1"
            + "\nOther versions:https://doi.org/10.6084/m9.figshare.6815705",
        ],
        # https://doi.org/10.6084/m9.figshare.6815699
        "dft_3d": [
            "https://ndownloader.figshare.com/files/38521619",
            "jdft_3d-12-12-2022.json",
            "Obtaining 3D dataset 76k ...",
            "https://www.nature.com/articles/s41524-020-00440-1"
            + "\nOther versions:https://doi.org/10.6084/m9.figshare.6815699",
        ],
        # https://doi.org/10.6084/m9.figshare.6815705
        "dft_2d_2021": [
            "https://ndownloader.figshare.com/files/26808917",
            "d2-3-12-2021.json",
            "Obtaining 2D dataset 1.1k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.6815699
        "dft_3d_2021": [
            "https://ndownloader.figshare.com/files/29204826",
            "jdft_3d-8-18-2021.json",
            "Obtaining 3D dataset 55k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.6815699
        "cfid_3d": [
            "https://ndownloader.figshare.com/files/29205201",
            "cfid_3d-8-18-2021.json",
            "Obtaining 3D dataset 55k ...",
            "https://www.nature.com/articles/s41524-020-00440-1"
            + "\nOther versions:https://doi.org/10.6084/m9.figshare.6815699",
        ],
        # https://doi.org/10.6084/m9.figshare.14213522
        "jff": [
            "https://ndownloader.figshare.com/files/28937793",
            # "https://ndownloader.figshare.com/files/26809760",
            "jff-7-24-2021.json",
            # "jff-3-12-2021.json",
            "Obtaining JARVIS-FF 2k ...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.21667874
        "alignn_ff_db": [
            "https://ndownloader.figshare.com/files/38522315",
            # "https://ndownloader.figshare.com/files/26809760",
            "id_prop.json",
            "Obtaining ALIGNN-FF training DB 300k ...",
            "https://doi.org/10.1039/D2DD00096B",
        ],
        "mp_3d_2020": [
            "https://ndownloader.figshare.com/files/26791259",
            "all_mp.json",
            "Obtaining Materials Project-3D CFID dataset 127k...",
            "https://doi.org/10.1063/1.4812323",
        ],
        # https://doi.org/10.6084/m9.figshare.14177630
        "megnet": [
            "https://ndownloader.figshare.com/files/26724977",
            "megnet.json",
            "Obtaining MEGNET-3D CFID dataset 69k...",
            "https://pubs.acs.org/doi/10.1021/acs.chemmater.9b01294",
        ],
        # https://doi.org/10.6084/m9.figshare.14745435
        "megnet2": [
            "https://ndownloader.figshare.com/files/28332741",
            "megnet-mp-2019-04-01.json",
            "Obtaining MEGNET-3D CFID dataset 133k...",
            "https://pubs.acs.org/doi/10.1021/acs.chemmater.9b01294",
        ],
        # https://doi.org/10.6084/m9.figshare.14745327
        "edos_pdos": [
            "https://ndownloader.figshare.com/files/29216859",
            "edos-up_pdos-elast_interp-8-18-2021.json",
            "Interpolated electronic total dos spin-up dataset 55k...",
            "https://www.nature.com/articles/s41524-020-00440-1",
        ],
        # https://doi.org/10.6084/m9.figshare.13054247
        "mp_3d": [
            "https://ndownloader.figshare.com/files/24979850",
            "CFID_mp_desc_data_84k.json",
            "Obtaining Materials Project-3D CFID dataset 84k...",
            "https://doi.org/10.1063/1.4812323",
        ],
        # https://doi.org/10.6084/m9.figshare.13055333
        "oqmd_3d": [
            "https://ndownloader.figshare.com/files/24981170",
            "CFID_OQMD_460k.json",
            "Obtaining OQMD-3D CFID dataset 460k...",
            "https://www.nature.com/articles/npjcompumats201510",
        ],
        # https://doi.org/10.6084/m9.figshare.14206169
        "oqmd_3d_no_cfid": [
            "https://ndownloader.figshare.com/files/26790182",
            "all_oqmd.json",
            "Obtaining OQMD-3D  dataset 800k...",
            "https://www.nature.com/articles/npjcompumats201510",
        ],
        # https://doi.org/10.6084/m9.figshare.14205083
        "twod_matpd": [
            "https://ndownloader.figshare.com/files/26789006",
            "twodmatpd.json",
            "Obtaining 2DMatPedia dataset 6k...",
            "https://www.nature.com/articles/s41597-019-0097-3",
        ],
        # https://doi.org/10.6084/m9.figshare.14213603
        "polymer_genome": [
            "https://ndownloader.figshare.com/files/26809907",
            "pgnome.json",
            "Obtaining Polymer genome 1k...",
            "https://www.nature.com/articles/sdata201612",
        ],
        "qm9_std_jctc": [
            "https://ndownloader.figshare.com/files/28715319",
            "qm9_std_jctc.json",
            "Obtaining QM9 standardized dataset 130k,"
            + "From https://doi.org/10.1021/acs.jctc.7b00577,+",
            "https://www.nature.com/articles/sdata201422",
        ],
        # https://doi.org/10.6084/m9.figshare.14827584
        # Use qm9_std_jctc instaed
        "qm9_dgl": [
            "https://ndownloader.figshare.com/files/28541196",
            "qm9_dgl.json",
            "Obtaining QM9 dataset 130k, from DGL...",
            "https://www.nature.com/articles/sdata201422",
        ],
        # https://doi.org/10.6084/m9.figshare.14912820.v1
        "cod": [
            "https://ndownloader.figshare.com/files/28715301",
            "cod_db.json",
            "Obtaining COD dataset 431k",
            "https://doi.org/10.1107/S1600576720016532",
        ],
        # Use qm9_std_jctc instaed
        "qm9": [
            "https://ndownloader.figshare.com/files/27627596",
            "qm9_data_cfid.json",
            "Obtaining QM9 dataset 134k...",
            "https://www.nature.com/articles/sdata201422",
        ],
        # https://doi.org/10.6084/m9.figshare.15127788
        "qe_tb": [
            "https://ndownloader.figshare.com/files/29070555",
            "jqe_tb_folder.json",
            "Obtaining QETB dataset 860k...",
            "https://arxiv.org/abs/2112.11585",
        ],
        # https://doi.org/10.6084/m9.figshare.14812050
        "omdb": [
            "https://ndownloader.figshare.com/files/28501761",
            "omdbv1.json",
            "Obtaining OMDB dataset 12.5k...",
            "https://doi.org/10.1002/qute.201900023",
        ],
        # https://doi.org/10.6084/m9.figshare.14812044
        "qmof": [
            "https://figshare.com/ndownloader/files/30972640",
            "qmof_db.json",
            "Obtaining QMOF dataset 20k...",
            "https://www.cell.com/matter/fulltext/S2590-2385(21)00070-9",
        ],
        # https://doi.org/10.6084/m9.figshare.15127758
        "hmof": [
            "https://figshare.com/ndownloader/files/30972655",
            "hmof_db_9_18_2021.json",
            "Obtaining hMOF dataset 137k...",
            "https://doi.org/10.1021/acs.jpcc.6b08729",
        ],
        # https://figshare.com/account/projects/100325/articles/14960157
        "c2db": [
            "https://ndownloader.figshare.com/files/28682010",
            "c2db_atoms.json",
            "Obtaining C2DB dataset 3.5k...",
            "https://iopscience.iop.org/article/10.1088/2053-1583/aacfc1",
        ],
        # https://doi.org/10.6084/m9.figshare.25256236
        "halide_peroskites": [
            "https://figshare.com/ndownloader/files/44619562",
            "halide_peroskites.json",
            "Obtaining halide perovskite dataset229...",
            "https://doi.org/10.1039/D1EE02971A",
        ],
        # https://figshare.com/account/projects/100325/articles/14962356
        "hopv": [
            "https://ndownloader.figshare.com/files/28814184",
            "hopv_15.json",
            "Obtaining HOPV15 dataset 4.5k...",
            "https://www.nature.com/articles/sdata201686",
        ],
        # https://figshare.com/account/projects/100325/articles/14962356
        "pdbbind_core": [
            "https://ndownloader.figshare.com/files/28874802",
            "pdbbind_2015_core.json",
            "Obtaining PDBBind dataset 195...",
            "https://doi.org/10.1093/bioinformatics/btu626",
        ],
        # https://doi.org/10.6084/m9.figshare.14812038
        "pdbbind": [
            "https://ndownloader.figshare.com/files/28816368",
            "pdbbind_2015.json",
            "Obtaining PDBBind dataset 11k...",
            "https://doi.org/10.1093/bioinformatics/btu626",
        ],
        # https://doi.org/10.6084/m9.figshare.21713885
        "snumat": [
            "https://ndownloader.figshare.com/files/38521736",
            "snumat.json",
            "Obtaining SNUMAT Hybrid functional dataset 10k...",
            "https://www.nature.com/articles/s41597-020-00723-8",
        ],
        # https://doi.org/10.6084/m9.figshare.13215308
        "aflow2": [
            "https://ndownloader.figshare.com/files/25453265",
            "CFID_AFLOW2.json",
            "Obtaining AFLOW-2 CFID dataset 400k...",
            "https://doi.org/10.1016/j.commatsci.2012.02.005",
        ],
        # https://doi.org/10.6084/m9.figshare.14211860
        "arXiv": [
            "https://ndownloader.figshare.com/files/26804795",
            "arXivdataset.json",
            "Obtaining arXiv dataset 1.8 million...",
            "https://www.kaggle.com/Cornell-University/arxiv",
        ],
        # https://doi.org/10.6084/m9.figshare.14211857
        "cord19": [
            "https://ndownloader.figshare.com/files/26804798",
            "cord19.json",
            "Obtaining CORD19 dataset 223k...",
            "https://github.com/usnistgov/cord19-cdcs-nist",
        ],
        # https://doi.org/10.6084/m9.figshare.22583677
        "ssub": [
            "https://figshare.com/ndownloader/files/40084921",
            "ssub.json",
            "Obtaining SSUB dataset 1726...",
            "https://github.com/wolverton-research-group/qmpy",
        ],
        # https://doi.org/10.6084/m9.figshare.22721047
        "mlearn": [
            # "https://figshare.com/ndownloader/files/40424156",
            "https://figshare.com/ndownloader/files/40357663",
            "mlearn.json",
            "Obtaining mlearn dataset 1730...",
            "https://github.com/materialsvirtuallab/mlearn",
        ],
        # https://doi.org/10.6084/m9.figshare.22814318
        "foundry_ml_exp_bandgaps": [
            "https://figshare.com/ndownloader/files/40557743",
            "foundry_ml_exp_bandgaps.json",
            "Obtaining foundry_ml_exp_bandgaps dataset 2069...",
            "https://foundry-ml.org/#/datasets/10.18126/wg3u-g8vu",
        ],
        # ToFix# https://doi.org/10.6084/m9.figshare.22815926
        # "mat_scholar_ner": [
        #    "https://figshare.com/ndownloader/files/40563593",
        #    "mat_scholar_ner.json",
        #    "Obtaining mat_scholar_ner dataset XYZ...",
        #    "https://pubs.acs.org/doi/10.1021/acs.jcim.9b00470",
        # ],
        # https://doi.org/10.6084/m9.figshare.22817633
        # Contains repeats
        "ocp10k": [
            "https://figshare.com/ndownloader/files/40566122",
            "ocp10k.json",
            "Obtaining OCP 10k train dataset, 59886...",
            "https://github.com/Open-Catalyst-Project/ocp",
        ],
        # https://doi.org/10.6084/m9.figshare.22817651
        "arxiv_summary": [
            "https://figshare.com/ndownloader/files/40566137",
            "arxiv_summary.json",
            "Obtaining arxiv summary cond.mat dataset 137927...",
            "https://github.com/usnistgov/chemnlp",
        ],
        # TODO:PubChem
        # https://doi.org/10.6084/m9.figshare.22975787
        "supercon_chem": [
            "https://figshare.com/ndownloader/files/40719260",
            "supercon_chem.json",
            "Obtaining supercon chem dataset 16414...",
            "https://www.nature.com/articles/s41524-018-0085-8",
        ],
        # https://doi.org/10.6084/m9.figshare.22976285
        "mag2d_chem": [
            "https://figshare.com/ndownloader/files/40720004",
            "mag2d_chem.json",
            "Obtaining magnetic 2D chem dataset 226...",
            "https://doi.org/10.24435/materialscloud:2019.0020/v1",
        ],
        # https://doi.org/10.6084/m9.figshare.23000573
        "vacancydb": [
            "https://figshare.com/ndownloader/files/40750811",
            "vacancydb.json",
            "Obtaining vacancy dataset 464...",
            "https://doi.org/10.1063/5.0135382",
        ],
        # https://doi.org/10.6084/m9.figshare.25832614
        "surfacedb": [
            "https://figshare.com/ndownloader/files/46355689",
            "surface_db_dd.json",
            "Obtaining vacancy dataset 607...",
            "https://doi.org/10.1039/D4DD00031E",
        ],
        # https://doi.org/10.6084/m9.figshare.25832614
        "interfacedb": [
            "https://figshare.com/ndownloader/files/46355692",
            "interface_db_dd.json",
            "Obtaining vacancy dataset 607...",
            "https://doi.org/10.1039/D4DD00031E",
        ],
        # Contains repeats
        # https://doi.org/10.6084/m9.figshare.23206193
        "ocp100k": [
            "https://figshare.com/ndownloader/files/40902845",
            "ocp100k.json",
            "Obtaining OCP100k dataset 149886...",
            "https://github.com/Open-Catalyst-Project/ocp",
        ],
        # https://doi.org/10.6084/m9.figshare.23250629
        "ocp_all": [
            "https://figshare.com/ndownloader/files/40974599",
            "ocp_all.json",
            "Obtaining OCPall dataset 510214...",
            "https://github.com/Open-Catalyst-Project/ocp",
        ],
        # https://doi.org/10.6084/m9.figshare.23225687
        "tinnet_N": [
            "https://figshare.com/ndownloader/files/40934285",
            "tinnet_N.json",
            "Obtaining TinNet Nitrogen dataset 329...",
            "https://github.com/hlxin/tinnet",
        ],
        # https://doi.org/10.6084/m9.figshare.23254151
        "tinnet_O": [
            "https://figshare.com/ndownloader/files/40978943",
            "tinnet_O.json",
            "Obtaining TinNet Oxygen dataset 747...",
            "https://github.com/hlxin/tinnet",
        ],
        # https://doi.org/10.6084/m9.figshare.23254154
        "tinnet_OH": [
            "https://figshare.com/ndownloader/files/40978949",
            "tinnet_OH.json",
            "Obtaining TinNet OH dataset 748...",
            "https://github.com/hlxin/tinnet",
        ],
        # https://doi.org/10.6084/m9.figshare.23909478
        "AGRA_O": [
            "https://figshare.com/ndownloader/files/41923284",
            "AGRA_O.json",
            "Obtaining AGRA Oxygen dataset 1000...",
            "https://github.com/Feugmo-Group/AGRA",
        ],
        # https://doi.org/10.6084/m9.figshare.23909478
        "AGRA_OH": [
            "https://figshare.com/ndownloader/files/41923287",
            "AGRA_OH.json",
            "Obtaining AGRA OH dataset 875...",
            "https://github.com/Feugmo-Group/AGRA",
        ],
        # https://doi.org/10.6084/m9.figshare.23909478
        "AGRA_CO": [
            "https://figshare.com/ndownloader/files/41923278",
            "AGRA_CO.json",
            "Obtaining AGRA CO dataset 193...",
            "https://github.com/Feugmo-Group/AGRA",
        ],
        # https://doi.org/10.6084/m9.figshare.23909478
        "AGRA_CHO": [
            "https://figshare.com/ndownloader/files/41923275",
            "AGRA_CHO.json",
            "Obtaining AGRA Oxygen dataset 214...",
            "https://github.com/Feugmo-Group/AGRA",
        ],
        # https://doi.org/10.6084/m9.figshare.23909478
        "AGRA_COOH": [
            "https://figshare.com/ndownloader/files/41923281",
            "AGRA_COOH.json",
            "Obtaining AGRA COOH dataset 280...",
            "https://github.com/Feugmo-Group/AGRA",
        ],
        # https://doi.org/10.6084/m9.figshare.21370572
        "supercon_3d": [
            "https://figshare.com/ndownloader/files/38307921",
            "jarvis_epc_data_figshare_1058.json",
            "Obtaining supercond. Tc dataset 1058...",
            "https://www.nature.com/articles/s41524-022-00933-1",
        ],
        # https://doi.org/10.6084/m9.figshare.21370572
        "supercon_2d": [
            "https://figshare.com/ndownloader/files/38950433",
            "jarvis_epc_data_2d.json",
            "Obtaining supercond. Tc dataset 161...",
            "https://doi.org/10.1021/acs.nanolett.2c04420",
        ],
        # https://doi.org/10.6084/m9.figshare.23267852
        "m3gnet_mpf": [
            "https://figshare.com/ndownloader/files/41009036",
            "m3gnet_mpf.json",
            "Obtaining m3gnet_mpf dataset 168917...",
            "https://github.com/materialsvirtuallab/m3gnet",
        ],
        # https://doi.org/10.6084/m9.figshare.23267852
        "m3gnet_mpf_1.5mil": [
            "https://figshare.com/ndownloader/files/47281519",
            "id_prop.json",
            "Obtaining m3gnet_mpf dataset 1.5mil...",
            "https://github.com/materialsvirtuallab/m3gnet",
        ],
        # https://doi.org/10.6084/m9.figshare.23531523
        "mxene275": [
            "https://figshare.com/ndownloader/files/41266233",
            "mxene275.json",
            "Obtaining mxene dataset 275...",
            "https://cmr.fysik.dtu.dk/c2db/c2db.html",
        ],
        # https://doi.org/10.6084/m9.figshare.26117998
        "cccbdb": [
            "https://figshare.com/ndownloader/files/47283808",
            "cccbdb.json",
            "Obtaining CCCBDB dataset 1333...",
            "https://cccbdb.nist.gov/",
        ],
        # https://doi.org/10.6084/m9.figshare.27174897
        "alex_pbe_hull": [
            "https://figshare.com/ndownloader/files/49622718",
            "alexandria_convex_hull_pbe_2023.12.29_jarvis_tools.json",
            "Obtaining Alexandria_DB PBE on hull 116k...",
            "https://alexandria.icams.rub.de/",
        ],
        # https://doi.org/10.6084/m9.figshare.27174897
        "alex_pbe_3d_all": [
            "https://figshare.com/ndownloader/files/49622946",
            "alexandria_pbe_3d_2024.10.1_jarvis_tools.json",
            "Obtaining Alexandria_DB PBE 3D all 5 million, large file...",
            "https://alexandria.icams.rub.de/",
        ],
        # https://doi.org/10.6084/m9.figshare.27174897
        "alex_pbe_2d_all": [
            "https://figshare.com/ndownloader/files/49622988",
            "alexandria_pbe_2d_2024.10.1_jarvis_tools.json",
            "Obtaining Alexandria_DB PBE 2D all 200k...",
            "https://alexandria.icams.rub.de/",
        ],
        # https://doi.org/10.6084/m9.figshare.27174897
        "alex_pbe_1d_all": [
            "https://figshare.com/ndownloader/files/49622991",
            "alexandria_pbe_1d_2024.10.1_jarvis_tools.json",
            "Obtaining Alexandria_DB PBE 1D all 100k...",
            "https://alexandria.icams.rub.de/",
        ],
        # https://doi.org/10.6084/m9.figshare.27174897
        "alex_scan_3d_all": [
            "https://figshare.com/ndownloader/files/49623090",
            "alexandria_scan_3d_2024.10.1_jarvis_tools.json",
            "Obtaining Alexandria_DB SCAN 3D all 500k...",
            "https://alexandria.icams.rub.de/",
        ],
        # https://doi.org/10.6084/m9.figshare.27174897
        "alex_pbesol_3d_all": [
            "https://figshare.com/ndownloader/files/49623096",
            "alexandria_ps_3d_2024.10.1_jarvis_tools.json",
            "Obtaining Alexandria_DB PBEsol 3D all 500k...",
            "https://alexandria.icams.rub.de/",
        ],
        # https://doi.org/10.6084/m9.figshare.13154159
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
    store_dir=None,
):
    """Get data with progress bar."""
    zfile = js_tag + ".zip"
    if store_dir is None:
        path = str(os.path.join(os.path.dirname(__file__), zfile))
    else:
        path = str(os.path.join(store_dir, zfile))

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


def data(dataset="dft_2d", store_dir=None):
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
    dat = get_request_data(js_tag=js_tag, url=url, store_dir=store_dir)
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
    from jarvis.analysis.stm.tersoff_hamann import TersoffHamannSTM

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
    from jarvis.io.wannier.outputs import WannierHam
    from jarvis.io.vasp.inputs import Poscar

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
    from jarvis.io.vasp.outputs import Vasprun
    from jarvis.io.wannier.outputs import WannierHam

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
    temp = np.zeros(w.nwan**2, dtype=complex)
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
