{ nixpkgs ? import ./nix/nixpkgs_version.nix }:
let
  pypkgs = nixpkgs.python36Packages;
  custodian = import ./nix/custodian.nix { inherit pypkgs; };
  fortranformat = import ./nix/fortranformat.nix { inherit pypkgs; };
  phonopy = import ./nix/phonopy.nix { inherit pypkgs; };
  pymatgen = pypkgs.pymatgen.overrideDerivation ( oldAttrs: rec {
    version = "2018.12.12";
    pname = "pymatgen";
    src = pypkgs.fetchPypi {
      inherit pname version;
      sha256 = "1isgwqxp24rd5i1z8w8cfj48bcx068g0h3jrk75n5p95dwrdy32j";
    };
  });
  ase = pypkgs.ase.overrideDerivation ( oldAttrs: rec {
    version = "3.11.0";
    pname = "ase";
    src = pypkgs.fetchPypi {
      inherit pname version;
      sha256 = "0kcfsx6wx2rdndy3463yngs4ap8zl0rjq6piwk24jnyrh6awqsyn";
    };
  });
in
  pypkgs.buildPythonPackage rec {
     pname = "jarvis";
     version = "dev";
     env = nixpkgs.buildEnv { name=pname; paths=buildInputs; };
     buildInputs = [
       pypkgs.pip
       pypkgs.python
       pypkgs.numpy
       pypkgs.scipy
       nixpkgs.pkgs.git
       pypkgs.matplotlib
       pypkgs.tkinter
       pymatgen
       ase
       pypkgs.pybtex
       custodian
       fortranformat
       pypkgs.networkx
       pypkgs.scikitlearn
       pypkgs.pytest
       pypkgs.h5py
       pypkgs.interruptingcow
       pypkgs.pybtex
       pypkgs.black
       pypkgs.pandas
     ];
     src=./.;
     doCheck=false;
     meta = {
       homepage = "https://github.com/usnistgov/jarvis";
       description = ''Joint Automated Repository for Various
       Integrated Simulations (JARVIS) is an integrated framework for
       computational science using density functional theory,
       classical force-field/molecular dynamics and
       machine-learning.'';
       version = version;
       license = nixpkgs.stdenv.lib.licenses.free;
     };
     catchConflicts=false;
     postShellHook = ''
       SOURCE_DATE_EPOCH=$(date +%s)
       export PYTHONUSERBASE=$PWD/.local
       export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
       export PYTHONPATH=$PYTHONPATH:$USER_SITE
       export PATH=$PATH:$PYTHONUSERBASE/bin

       ## To build a python package from pypi use
       # pip install --user package_name
     '';
  }
