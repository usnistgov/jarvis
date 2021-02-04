"""
Module for obtaining webpage data.

All the webpages are formatted in XML format.
They are converted into python dict from which
appropriate keys and values can be obtained.
Use 'list_keys' functions to find types
of data stored in each XML document.
"""
import numpy as np
import os
import requests
import tempfile
from jarvis.core.utils import xml_to_dict


class Webpage(object):
    """Get content of JARVIS-DFT and JARVIS-FF webpages."""

    def __init__(self, data=[], jid=""):
        """Initialize the class."""
        # JVASP- for JARVIS-DFT IDs
        # JLMP- for JARVIS-FF IDs

        self.jid = jid

        if "JVASP-" in jid:
            url = (
                "https://www.ctcms.nist.gov/~knc6/static/JARVIS-DFT/"
                + jid
                + ".xml"
            )
        if "JLMP-" in jid:
            url = (
                "https://www.ctcms.nist.gov/~knc6/static/JARVIS-FF/"
                + jid
                + ".xml"
            )
        self.url = url
        self.data = data
        if self.data == []:
            data = self.to_dict()
            self.data = data

    def to_dict(self):
        """Get dictionary for webpages."""
        dat = requests.get(self.url).text
        fd, path = tempfile.mkstemp()
        with os.fdopen(fd, "w") as tmp:
            tmp.write(dat)
        data = xml_to_dict(path)
        return data

    def list_keys(self):
        """Get all main keys from webpages."""
        print("keys=", self.data["basic_info"].keys())
        return self.data["basic_info"].keys()

    def get_dft_electron_dos(self):
        """Get JARVIS-DFT MAIN-RELAX DOS."""
        return self.data["basic_info"]["main_relax_info"]["main_relax_dos"]

    def get_dft_phonon_dos(self):
        """Get JARVIS-DFT MAIN-ELAST finite-diff DOS."""
        return self.data["basic_info"]["main_elastic"]["main_elastic_info"]

    def get_dft_mbj_dielectric_function(self):
        """Get JARVIS-DFT TBmBJ dielectric function."""
        info = {}
        energies = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "energies"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_xx = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "real_1"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_yy = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "real_2"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_zz = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "real_3"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_xy = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "real_4"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_yz = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "real_5"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_zx = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "real_6"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_xx = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "imag_1"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_yy = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "imag_2"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_zz = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "imag_3"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_xy = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "imag_4"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_yz = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "imag_5"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_zx = np.array(
            self.data["basic_info"]["main_optics_mbj"]["main_optics_mbj_info"][
                "imag_6"
            ]
            .strip("'")
            .split(","),
            dtype="float",
        )
        info["energies"] = energies
        info["real_xx"] = real_part_xx
        info["real_yy"] = real_part_yy
        info["real_zz"] = real_part_zz
        info["real_xy"] = real_part_xy
        info["real_yz"] = real_part_yz
        info["real_zx"] = real_part_zx
        info["imag_xx"] = imag_part_xx
        info["imag_yy"] = imag_part_yy
        info["imag_zz"] = imag_part_zz
        info["imag_xy"] = imag_part_xy
        info["imag_yz"] = imag_part_yz
        info["imag_zx"] = imag_part_zx
        return info

    def get_dft_semilocal_dielectric_function(self):
        """Get JARVIS-DFT semilocal dielectric function."""
        # OptB88vdW/LDA/PBE

        info = {}
        energies = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["energies"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_xx = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["real_1"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_yy = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["real_2"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_zz = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["real_3"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_xy = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["real_4"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_yz = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["real_5"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        real_part_zx = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["real_6"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_xx = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["imag_1"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_yy = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["imag_2"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_zz = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["imag_3"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_xy = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["imag_4"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_yz = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["imag_5"]
            .strip("'")
            .split(","),
            dtype="float",
        )
        imag_part_zx = np.array(
            self.data["basic_info"]["main_optics_semilocal"][
                "main_optics_info"
            ]["imag_6"]
            .strip("'")
            .split(","),
            dtype="float",
        )

        info["energies"] = energies
        info["real_xx"] = real_part_xx
        info["real_yy"] = real_part_yy
        info["real_zz"] = real_part_zz
        info["real_xy"] = real_part_xy
        info["real_yz"] = real_part_yz
        info["real_zx"] = real_part_zx
        info["imag_xx"] = imag_part_xx
        info["imag_yy"] = imag_part_yy
        info["imag_zz"] = imag_part_zz
        info["imag_xy"] = imag_part_xy
        info["imag_yz"] = imag_part_yz
        info["imag_zx"] = imag_part_zx
        return info


"""
if __name__ == "__main__":
    w = Webpage(jid="JVASP-1002")
    info_mbj = w.get_dft_mbj_dielectric_function()
    # print(info_mbj)
    w.list_keys()
"""
