{ pypkgs }:
pypkgs.buildPythonPackage rec {
  version = "0.2.5";
  pname = "fortranformat";
  src = pypkgs.fetchPypi {
    inherit pname version;
    sha256 = "1yfkki9yimznd7724jw3y9yz0msdzcsg709s7ia70ylw28gvqpvb";
  };
  doCheck=false;
  buildInputs = [];
}
