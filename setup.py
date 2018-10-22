import glob
import os

from setuptools import setup, find_packages

JARVIS_DIR = os.path.dirname(os.path.abspath(__file__))

setup(
    name="jarvis-tools",
    version="2018.10.22.1",
    install_requires=["custodian>=1.0.4","ase==3.11.0","pymatgen==2017.8.4",
                 "pybtex==0.21","phonopy==1.10.8","fortranformat==0.2.5",
                 "pandas==0.23.4","scikit-learn==0.18.1","numpy==1.15.2","scipy==1.1.0","lightgbm==2.1.0",]

    install_requires=["custodian>=1.0.4","ase==3.11.0","pymatgen==2017.8.4","sphinx==1.8.1","pybtex==0.21","phonopy==1.10.8","fortranformat==0.2.5","pandas==0.23.4","scikit-learn==0.18.1","lightgbm==2.1.0"],
    extras_require={"babel": ["openbabel", "pybel"],
                    "remote": ["fabric"],
                    "doc": ["sphinx>=1.3.1", "sphinx-rtd-theme>=0.1.8"]
                    },
    author="Kamal Choudhary",
    author_email = "kamal.choudhary@nist.gov",
    description=(
        "High throughput computation with density functional theory, molecular dynamics and machine learning"),
    license="MIT",
    url="https://github.com/usnistgov/jarvis",
    packages=find_packages(),
    long_description=open(
        os.path.join(os.path.dirname(__file__), 'README.rst')).read(),
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
    #scripts=glob.glob(os.path.join(JARVIS_DIR,  "*"))
)

