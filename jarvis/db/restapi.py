"""
To access data in JARVIS-API, you do not need to request an account
However, to upload your data, please request an account at https://jarvis.nist.gov/
Then install MDCS-api-tools by:

git clone https://github.com/knc6/MDCS-api-tools.git
cd MDCS-api-tools
python setup.py develop

You'll require a schema ID and a user-name to go through examples below
"""

from jarvis.analysis.structure.spacegroup import Spacegroup3D
from jarvis.db.figshare import data
import os, glob, json, sys

try:
    import mdcs
    from mdcs import curate, explore
except:
    print("Please install MDCS")
    pass
import numpy as np
from jarvis.io.vasp.outputs import Vasprun
from jarvis.core.atoms import Atoms
from xml.etree.ElementTree import Element, SubElement, Comment, tostring, ElementTree

try:
  f = open("passwrd.txt", "r")  # path to your passowrd file
  passd = f.read().splitlines()[0]
  f.close()
except:
   print ('Passowrd is not provided')
   pass

user = "abc"  # your username
schema = "5a7f2a872887d500ab0c0d02"


def data_json(
    energy="na",
    typ="na",
    formula="na",
    sgp="na",
    name="na",
    ref="na",
    func="na",
    elem="na",
    encut="na",
    kpoints="na",
    el_tens="na",
    KV="na",
    GV="na",
    m_eg="na",
    b_eg="na",
    op_eg="na",
    mbj_eg="na",
    en_arr="na",
    realx_arr="na",
    imagx_arr="na",
    realy_arr="na",
    imagy_arr="na",
    realz_arr="na",
    imagz_arr="na",
    men_arr="na",
    mrealx_arr="na",
    mimagx_arr="na",
    mrealy_arr="na",
    mimagy_arr="na",
    mrealz_arr="na",
    mimagz_arr="na",
    struct="na",
    other="na",
    curate_xml=False
):
    """
    This is an example of how to upload JARVIS data to the API. 
    It converts input data to XML based on a template schema which must be available in the JARVIS-API.
    After the conversion, the XML could be uploaded to the API.
    Any user with registered username and password can upload.
    If they prefer to upload using their own schema, contact the developers.
    See an example: https://jarvis.nist.gov/data?id=5df7f05ceaf3b300338bf83f
    Args:

         metadata with a calculation
     
         curate_xml: for debugging set it to False, and produce XML files
   
    """

    top = Element("JARVIS-DFT")
    child = SubElement(top, "energy-ev")
    child.text = str(energy)
    child = SubElement(top, "formula")
    child.text = str(formula)
    child = SubElement(top, "space-group")
    child.text = str(sgp)
    child = SubElement(top, "JVID")
    child.text = str(name)
    child = SubElement(top, "reference")
    child.text = str(ref)
    child = SubElement(top, "calc_type")
    child.text = str(typ)
    child = SubElement(top, "functional")
    child.text = str(func)
    child = SubElement(top, "elements")
    child.text = str(elem)
    child = SubElement(top, "encut-ev")
    child.text = str(encut)
    child = SubElement(top, "kpoints")
    child.text = str(kpoints)
    child = SubElement(top, "elastic_tensor-gpa")
    child.text = str(el_tens)
    child = SubElement(top, "kv-gpa")
    child.text = str(KV)
    child = SubElement(top, "gv-gpa")
    child.text = str(GV)
    child = SubElement(top, "scf_eg-ev")
    child.text = str(m_eg)
    child = SubElement(top, "brill_eg-ev")
    child.text = str(b_eg)
    child = SubElement(top, "optics_eg-ev")
    child.text = str(op_eg)
    child = SubElement(top, "mbj-eg_ev")
    child.text = str(mbj_eg)
    child = SubElement(top, "opt_energy_arr-ev")
    child.text = str(en_arr)
    child = SubElement(top, "realx_arr")
    child.text = str(realx_arr)
    child = SubElement(top, "imagx_arr")
    child.text = str(imagx_arr)
    child = SubElement(top, "realy_arr")
    child.text = str(realy_arr)
    child = SubElement(top, "imagy_arr")
    child.text = str(imagy_arr)
    child = SubElement(top, "realz_arr")
    child.text = str(realz_arr)
    child = SubElement(top, "imagz_arr")
    child.text = str(imagz_arr)
    child = SubElement(top, "mbj_energy_arr-ev")
    child.text = str(men_arr)
    child = SubElement(top, "mbj_realx_arr")
    child.text = str(mrealx_arr)
    child = SubElement(top, "mbj_imagx_arr")
    child.text = str(mimagx_arr)
    child = SubElement(top, "mbj_realy_arr")
    child.text = str(mrealy_arr)
    child = SubElement(top, "mbj_imagy_arr")
    child.text = str(mimagy_arr)
    child = SubElement(top, "mbj_realz_arr")
    child.text = str(mrealz_arr)
    child = SubElement(top, "mbj_imagz_arr")
    child.text = str(mimagz_arr)
    child = SubElement(top, "structure")
    child.text = str(struct)
    child = SubElement(top, "other")
    child.text = str(other)
    filename = str(name) + str(".xml")
    ElementTree(top).write(filename)
    if curate_xml==True:
      curate(
        filename, filename, schema, "https://jarvis.nist.gov/", user, passd, cert=False
      )


def get_record(file="JVASP-1002.xml"):
   """
   This is an example of how to get a particular jarvis-id document, say JVASP-1002.xml

   Args:

       file: name of the file

   Returns:

         data in json format
   """

   r = explore.select(
        "https://jarvis.nist.gov/", user, passd, cert=False, title=file, format="json"
   )
   return r


def delete_all(file=""):
    """
    Caution, this deletes your documents
    """
    r = explore.select_all(
        "https://jarvis.nist.gov/", user, passd, cert=False, format="json"
    )
    for i in r:
        id = i["_id"]
        explore.delete(id, "https://jarvis.nist.gov/", user, passd, cert=False)


def upload_sample_data(curate_xml=False):
    """
    Generate and upload XML files
    Set curate_xml==True to upload all the XML documents
    """
    d = data('dft_2d')

    count = 0
    for i in d[0:1]:
        filname = str(i["jid"]) + str(".xml")
        if not os.path.exists(filname):
            count = count + 1
            energy = str(i["optb88vdw_total_energy"]) + str(",") + str(i["formation_energy_peratom"])
            atoms = Atoms.from_dict(i["atoms"])
            formula = str(atoms.composition.reduced_formula)
            sgp = str(Spacegroup3D(atoms).space_group_symbol)
            name = str(i["jid"])
            print(name)
            ref = str("")
            func = str("OptB88vdW")
            elem = ""
            species = atoms.elements
            for j in species:
                elem = str(elem) + str(j) + str("-")
            encut = str(i["encut"])
            kpoints = (
                str(i["kpoints_array"][0])
                + str("x")
                + str(i["kpoints_array"][1])
                + str("x")
                + str(i["kpoints_array"][2])
            )
            el_tens = 'na'
            try:
                el_tens = str(','.join(map(str,np.array(i["elastic_tensor"]).flatten())))
            except:
                pass
            KV = str(i["bulk_modulus_kv"])
            GV = str(i["shear_modulus_gv"])
            op_eg = str(i["optb88vdw_bandgap"])
            mbj_eg = str(i["mbj_bandgap"])
            realx_arr = str(i["epsx"])
            mrealx_arr = str(i["mepsx"])
            realy_arr = str(i["epsy"])
            mrealy_arr = str(i["mepsy"])
            realz_arr = str(i["epsz"])
            mrealz_arr = str(i["mepsz"])
            typ = str("2D") #3D
            other = str(
                "Citation: 1) DOI:10.1038/s41598-017-05402-0, 2) DOI: 10.1038/sdata.2018.82, 3) arXiv:1804.01033v2 "
            )
            struct = atoms.get_string()
            data_json(
                other=other,
                energy=energy,
                typ=typ,
                formula=formula,
                sgp=sgp,
                name=name,
                func=func,
                elem=elem,
                encut=encut,
                kpoints=kpoints,
                el_tens=el_tens,
                KV=KV,
                GV=GV,
                op_eg=op_eg,
                mbj_eg=mbj_eg,
                realx_arr=realx_arr,
                mrealx_arr=mrealx_arr,
                realy_arr=realy_arr,
                mrealy_arr=mrealy_arr,
                realz_arr=realz_arr,
                mrealz_arr=mrealz_arr,
                struct=struct,
                curate_xml=curate_xml
            )
"""
if __name__ == "__main__":
     upload_sample_data(curate_xml = False)
     #x=get_record(file='JVASP-48137.xml')[0]['content']['JARVIS-DFT']['structure']
     #print (x)
"""
