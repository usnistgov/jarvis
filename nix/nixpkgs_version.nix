let
  inherit (import <nixpkgs> {}) fetchFromGitHub;
  nixpkgs_download = fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "6a3f5bcb061e1822f50e299f5616a0731636e4e7"; #"7f35ed9df40f12a79a242e6ea79b8a472cf74d42";
    sha256 = "1ib96has10v5nr6bzf7v8kw7yzww8zanxgw2qi1ll1sbv6kj6zpd"; #"1wr6dzy99rfx8s399zjjjcffppsbarxl2960wgb0xjzr7v65pikz";
    # rev = "dcf40f7c24eec1160e6433b6644d3e2dd268e417";
    # sha256 = "0ajwr7n72njkj87asr29bs94pyg61ix6fx7a239n9ik6bisg639y";
  };
in
  import nixpkgs_download {}
