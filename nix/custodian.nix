{ pypkgs }:
pypkgs.buildPythonPackage rec {
  version = "2018.8.10";
  pname = "custodian";
  src = pypkgs.fetchPypi {
    inherit pname version;
    sha256 = "0s91lf43kx3450lqyw6adb8mp5r3mj9wkr450jvpxvxvc51zxp9d";
  };
  doCheck=false;
  buildInputs = [
    pypkgs.monty
    pypkgs.pyyaml
  ];
}
