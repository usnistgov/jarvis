"""Provide LAMMPS module .mod files."""

import os
import shutil


class GenericInputs(object):
    """Provide generic LAMMPS inputs such as relxation tolerance etc."""

    def __init__(self, input_name=""):
        """Intialize with filename."""
        self.input_name = input_name

    def elastic_general(self, path="."):
        """Provide filename for elastic constant calculations."""
        tmp_file = str(os.path.join(os.path.dirname(__file__), "inelast.mod"))
        shutil.copy2(tmp_file, path)

    def relax(self, path="."):
        """Provide filename for relaxation calculations."""
        tmp_file = str(os.path.join(os.path.dirname(__file__), "relax.mod"))
        shutil.copy2(tmp_file, path)

    def run0(self, path="."):
        """Provide filename for one-step calculations."""
        tmp_file = str(os.path.join(os.path.dirname(__file__), "run0.mod"))
        shutil.copy2(tmp_file, path)
