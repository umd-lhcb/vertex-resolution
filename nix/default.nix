{
  stdenv
, root
, cxxopts
, libyamlcpp
, boost
}:

stdenv.mkDerivation {
  pname = "vertex-resolution";
  version = "0.1.0";

  src = builtins.path { path = ./..; name = "vertex-resolution"; };

  buildInputs = [ root cxxopts boost ];

  buildPhase = "make";

  installPhase = ''
    mkdir -p $out/bin
    cp bin/* $out/bin
  '';
}
