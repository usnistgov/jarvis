{ pypkgs }:
pypkgs.buildPythonPackage rec {
  version = "1.13.2.107";
  pname = "phonopy";
  src = pypkgs.fetchPypi {
    inherit pname version;
    sha256 = "01grl0h2c7lzqbxzrs7lrzxjk402bnpa4vrdvfrpm6lbr4l6gw3j";
  };
  doCheck = false;
  buildInputs = [
    pypkgs.numpy
    pypkgs.matplotlib
    pypkgs.tkinter
    pypkgs.h5py
    pypkgs.pyyaml
  ];
}
