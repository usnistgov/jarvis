"""Module to prepare input files."""


def get_Zeff():
    """Get valence electrons."""
    return dict(
        Ag=19,
        Al=3,
        Ar=8,
        As=5,
        Au=19,
        B=3,
        Be=2,
        Bi=5,
        Br=7,
        C=4,
        Ca=10,
        Cl=7,
        Co=17,
        Cr=14,
        Cu=19,
        F=7,
        Fe=16,
        Ga=3,
        Ge=4,
        H=1,
        He=2,
        I=7,
        Ir=17,
        K=9,
        Kr=8,
        Li=1,
        Mg=2,
        Mn=15,
        Mo=14,
        N=5,
        Na=1,
        Ne=8,
        Ni=18,
        O=6,
        P=5,
        Pd=18,
        S=6,
        Sc=11,
        Se=6,
        Si=4,
        Tb=19,
        Te=6,
        Ti=12,
        V=13,
        W=14,
        Zn=20,
    )


def get_pseudo_dft_dict():
    """Get PSPs for DFT."""
    return {
        "Ag": "Ag.ccECP.AREP.upf",  # Ag.ccECP.SOREP.upf
        "Al": "Al.ccECP.upf",
        "Ar": "Ar.ccECP.upf",
        "As": "As.ccECP.upf",
        "Au": "Au.ccECP.AREP.upf",  # Au.ccECP.SOREP.upf
        "B": "B.ccECP.upf",
        "Be": "Be.ccECP.upf",
        "Bi": "Bi.ccECP.AREP.upf",  # Bi.ccECP.SOREP.upf
        "Br": "Br.ccECP.upf",
        "C": "C.ccECP.upf",
        "Ca": "Ca.ccECP.upf",
        "Cl": "Cl.ccECP.upf",
        "Co": "Co.ccECP-soft.upf",  # Co.opt.upf
        "Cr": "Cr.ccECP.upf",  # Cr.opt.upf
        "Cu": "Cu.ccECP-soft.upf",  # Cu.opt.upf
        "F": "F.ccECP.upf",
        "Fe": "Fe.ccECP-soft.upf",  # Fe.opt.upf
        "Ga": "Ga.ccECP.upf",
        "Ge": "Ge.ccECP.upf",
        "H": "H.ccECP.upf",
        "He": "He.ccECP.upf",
        "I": "I.ccECP.AREP.upf",  # I.ccECP.SOREP.upf
        "Ir": "Ir.ccECP.AREP.upf",  # Ir.ccECP.SOREP.upf
        "K": "K.ccECP.upf",
        "Kr": "Kr.ccECP.upf",
        "Li": "Li.ccECP.upf",
        "Mg": "Mg.ccECP.upf",
        "Mn": "Mn.ccECP.upf",  # Mn.opt.upf
        "Mo": "Mo.ccECP.AREP.upf",  # Mo.ccECP.SOREP.upf
        "N": "N.ccECP.upf",
        "Na": "Na.ccECP.upf",
        "Ne": "Ne.ccECP.upf",
        "Ni": "Ni.ccECP-soft.upf",  # Ni.opt.upf
        "O": "O.ccECP.upf",
        "P": "P.ccECP.upf",
        "Pd": "Pd.ccECP.AREP.upf",  # Pd.ccECP.SOREP.upf
        "S": "S.ccECP.upf",
        "Sc": "Sc.ccECP-soft.upf",  # Sc.opt.upf
        "Se": "Se.ccECP.upf",
        "Si": "Si.ccECP.upf",
        "Tb": "Tb.ccECP.AREP.upf",  # Tb.ccECP.SOREP.upf
        "Te": "Te.ccECP.AREP.upf",  # Te.ccECP.SOREP.upf
        "Ti": "Ti.ccECP-soft.upf",  # Ti.opt.upf
        "V": "V.ccECP-soft.upf",  # V.opt.upf
        "W": "W.ccECP.AREP.upf",  # W.ccECP.SOREP.upf
        "Zn": "Zn.ccECP-soft.upf",  # Zn.opt.upf
    }


def get_pseudo_qmc_dict():
    """Get PSPs for QMC."""
    return {
        "Ag": "Ag.ccECP.AREP.xml",  # Ag.ccECP.SOREP.xml
        "Al": "Al.ccECP.xml",
        "Ar": "Ar.ccECP.xml",
        "As": "As.ccECP.upf",
        "Au": "Au.ccECP.AREP.xml",  # Au.ccECP.SOREP.xml
        "B": "B.ccECP.xml",
        "Be": "Be.ccECP.xml",
        "Bi": "Bi.ccECP.AREP.xml",  # Bi.ccECP.SOREP.xml
        "Br": "Br.ccECP.xml",
        "C": "C.ccECP.xml",
        "Ca": "Ca.ccECP.xml",
        "Cl": "Cl.ccECP.xml",
        "Co": "Co.ccECP-soft.xml",  # Co.opt.xml
        "Cr": "Cr.ccECP.xml",  # Cr.opt.xml
        "Cu": "Cu.ccECP-soft.xml",  # Cu.opt.xml
        "F": "F.ccECP.xml",
        "Fe": "Fe.ccECP-soft.xml",  # Fe.opt.xml
        "Ga": "Ga.ccECP.xml",
        "Ge": "Ge.ccECP.xml",
        "H": "H.ccECP.xml",
        "He": "He.ccECP.xml",
        "I": "I.ccECP.AREP.xml",  # I.ccECP.SOREP.xml
        "Ir": "Ir.ccECP.AREP.xml",  # Ir.ccECP.SOREP.xml
        "K": "K.ccECP.xml",
        "Kr": "Kr.ccECP.xml",
        "Li": "Li.ccECP.xml",
        "Mg": "Mg.ccECP.xml",
        "Mn": "Mn.ccECP.xml",  # Mn.opt.xml
        "Mo": "Mo.ccECP.AREP.xml",  # Mo.ccECP.SOREP.xml
        "N": "N.ccECP.xml",
        "Na": "Na.ccECP.xml",
        "Ne": "Ne.ccECP.xml",
        "Ni": "Ni.ccECP-soft.xml",  # Ni.opt.xml
        "O": "O.ccECP.xml",
        "P": "P.ccECP.upf",
        "Pd": "Pd.ccECP.AREP.xml",  # Pd.ccECP.SOREP.xml
        "S": "S.ccECP.xml",
        "Sc": "Sc.ccECP-soft.xml",  # Sc.opt.xml
        "Se": "Se.ccECP.xml",
        "Si": "Si.ccECP.xml",
        "Tb": "Tb.ccECP.AREP.xml",  # Tb.ccECP.SOREP.xml
        "Te": "Te.ccECP.AREP.xml",  # Te.ccECP.SOREP.xml
        "Ti": "Ti.ccECP-soft.xml",  # Ti.opt.xml
        "V": "V.ccECP-soft.xml",  # V.opt.xml
        "W": "W.ccECP.AREP.xml",  # W.ccECP.SOREP.xml
        "Zn": "Zn.ccECP-soft.xml",  # Zn.opt.xml
    }
