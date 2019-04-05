{ nixpkgs ? import ./nixpkgs_version.nix }:
let
  pypkgs = nixpkgs.python36Packages;
in
  pypkgs.pymatgen.overrideDerivation ( oldAttrs: rec {
    version = "2018.12.12";
    pname = "pymatgen";
    src = pypkgs.fetchPypi {
      inherit pname version;
      sha256 = "1isgwqxp24rd5i1z8w8cfj48bcx068g0h3jrk75n5p95dwrdy32j";
    };
  })
