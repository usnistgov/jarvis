from jarvis.db.figshare import (
    data,
    get_ff_eneleast,
    get_wann_electron,
    get_wann_phonon,
    get_hk_tb,
    make_stm_from_prev_parchg,
    get_stm_2d_dataset,
)
from jarvis.db.webpages import Webpage


def test_figshare_download():

    data_2d = data(dataset="dft_2d")
    # data_3d = data(dataset="dft_3d")
    # data_ml = data(dataset="cfid_3d")
    data_ff = get_ff_eneleast()
    print("2d", len(data_2d))
    # print("3d", len(data_3d))
    # print("cfid3d", len(data_ml))
    print("ff", len(data_ff))
    # assert (len(data_2d), len(data_3d), len(data_ml), len(data_ff)) == (
    #    1070,
    #    36099,
    #    35984,
    #    3291,
    # )
    w, ef, atoms = get_wann_electron()
    hk = get_hk_tb(w=w)
    w, atoms = get_wann_phonon()
    w = Webpage(jid="JVASP-1002")
    info_mbj = w.get_dft_mbj_dielectric_function()
    info_opt = w.get_dft_semilocal_dielectric_function()
    info_pdos = w.get_dft_phonon_dos()
    info_edos = w.get_dft_electron_dos()
    k = w.list_keys()
    # dat = data(dataset="megnet")
    # dat = data(dataset="mp_3d")
    # dat = data(dataset="mp_3d_2020")
    # dat = data(dataset="qm9")
    # dat = data(dataset="aflow1")
    # dat = data(dataset="aflow2")
    # dat = data(dataset="oqmd_3d")
    # dat = data(dataset="twod_matpd")
    # dat = data(dataset="oqmd_3d_no_cfid")
    p, n = get_stm_2d_dataset()

    make_stm_from_prev_parchg()


# test_figshare_download()
