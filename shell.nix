{ pkgs ? (import (builtins.fetchTarball {
    url = "https://github.com/NixOS/nixpkgs/archive/20.03.tar.gz";
    sha256 = "0182ys095dfx02vl2a20j1hz92dx3mfgz2a6fhn31bqlp1wa8hlq";
  }) {}) }:
let
  pypkgs = pkgs.python3Packages;
in
  pypkgs.buildPythonPackage rec {
     pname = "jarvis-tools";
     version = "dev";
     nativeBuildInputs = with pypkgs; [
       joblib
       flask
       numpy
       phonopy
       scipy
       matplotlib
       spglib
       requests
       toolz
       pytest
       bokeh
       networkx
       xmltodict
     ];
     src=builtins.filterSource (path: type: type != "directory" || baseNameOf path != ".git") ./.;
     preShellHook = ''
       SOURCE_DATE_EPOCH=$(date +%s)
       export PYTHONUSERBASE=$PWD/.local
       export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
       export PYTHONPATH=$PYTHONPATH:$USER_SITE
       export PATH=$PATH:$PYTHONUSERBASE/bin

       ## To build a python package from pypi use
       pip install --user scikit-learn pandas lightgbm 
     '';
  }
