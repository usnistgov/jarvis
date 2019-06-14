{ nixpkgs, pypkgs }:
pypkgs.buildPythonPackage rec {
  name = "sfepy_${version}";
  version = "2018.2";
  src = nixpkgs.fetchurl {
    url="https://github.com/sfepy/sfepy/archive/release_${version}.tar.gz";
    sha256 = "16g51rnhfkfw1xs8qja837fbhqdnvfpnhy69jb1rvy5wd12qiqk9";
  };
  doCheck = false;
  buildInputs = [
    pypkgs.numpy
    pypkgs.cython
    pypkgs.scipy
    pypkgs.matplotlib
    pypkgs.pyparsing
    pypkgs.tables
  ];
  catchConflicts = false;
}
