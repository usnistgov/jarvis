import os
import shutil


class GenericInputs(object):
    def __init__(self, input_name=""):
        self.input_name = input_name

    def elastic_general(self, path="."):
        tmp_file = str(os.path.join(os.path.dirname(__file__), "inelast.mod"))
        shutil.copy2(tmp_file, path)

    def relax(self, path="."):
        tmp_file = str(os.path.join(os.path.dirname(__file__), "relax.mod"))
        shutil.copy2(tmp_file, path)

    def run0(self, path="."):
        tmp_file = str(os.path.join(os.path.dirname(__file__), "run0.mod"))
        shutil.copy2(tmp_file, path)
