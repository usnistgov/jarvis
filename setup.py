import glob
import os

from setuptools import setup, find_packages

JARVIS_DIR = os.path.dirname(os.path.abspath(__file__))

base_dir = os.path.dirname(__file__)
with open(os.path.join(base_dir, "README.rst")) as f:
    long_d = f.read()

setup(
    name="jarvis-tools",
    version="2019.07.30",

    long_description=long_d,
    long_description_content_type='text/markdown',
    install_requires=[
        "numpy==1.16.3",
        "scipy==1.2.1",
        "pymatgen>=2018.12.12",
        "custodian",
        "ase==3.11.0",
        "scikit-learn",
        "interruptingcow>=0.8",
        "pybtex>=0.21",
        "blob",
        "toolz",
        "fortranformat>=0.2.5",
      
    ],
    package_data={"jarvis.sklearn":['element_charge.json','Elements.json']},
    extras_require={
        "babel": ["openbabel", "pybel"],
        "doc": ["sphinx>=1.3.1", "sphinx-rtd-theme>=0.1.8"],
    },
    author="Kamal Choudhary",
    author_email="kamal.choudhary@nist.gov",
    description=(
        "High throughput computation with density functional theory, molecular dynamics and machine learning. https://jarvis.nist.gov/"
    ),
    license="MIT",
    url="https://github.com/usnistgov/jarvis",
    packages=find_packages(),
    #long_description=open(os.path.join(os.path.dirname(__file__), "README.rst")).read(),
    classifiers=[
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
    ],
    # scripts=glob.glob(os.path.join(JARVIS_DIR,  "*"))
)
