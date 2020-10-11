import glob
import os

from setuptools import setup, find_packages

JARVIS_DIR = os.path.dirname(os.path.abspath(__file__))

base_dir = os.path.dirname(__file__)
with open(os.path.join(base_dir, "README.rst")) as f:
    long_d = f.read()

setup(
    name="jarvis-tools",

    version="2020.10.10",

    long_description=long_d,
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "spglib",
        "joblib",
        "requests",
        "toolz",
        "xmltodict",
        "pytest",
    ],
    package_data={
        "jarvis.core": ["Elements.json", "element_charge.json"],
        "jarvis.tasks.lammps.templates": [
            "displace.mod",
            "inelastcomb.mod",
            "inelast_min.mod",
            "inelast.mod",
            "inelast_nobox.mod",
            "inelastreax.mod",
            "relax.mod",
            "run0.mod",
        ],
        "jarvis.io.vasp": ["default_potcars.json"],
        "jarvis.analysis.solarefficiency": ["am1.5G.dat"],
        "jarvis.io.wannier": ["default_semicore.json"],
    },
    extras_require={
        "ai": [
            "torch",
            "keras",
            "tensorflow",
            "scikit-learn",
            "flask",
            "pandas",
        ],
        "babel": ["openbabel", "pybel"],
        "doc": ["sphinx>=1.3.1", "sphinx-rtd-theme>=0.1.8"],
    },
    author="Kamal Choudhary",
    author_email="kamal.choudhary@nist.gov",
    description=(
        "jarvis-tools: an open-source software package for data-driven atomistic materials design. https://jarvis.nist.gov/"
    ),
    license="NIST",
    url="https://github.com/usnistgov/jarvis",
    packages=find_packages(),
    # long_description=open(os.path.join(os.path.dirname(__file__), "README.rst")).read(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    # scripts=glob.glob(os.path.join(JARVIS_DIR,  "*"))
)
