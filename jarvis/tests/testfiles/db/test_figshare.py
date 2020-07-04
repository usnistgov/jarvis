from jarvis.db.figshare import data, get_ff_eneleast


def test_figshare_download():

    data_2d = data(dataset="dft_2d")
    data_3d = data(dataset="dft_3d")
    data_ml = data(dataset="cfid_3d")
    data_ff = get_ff_eneleast()
    print("2d", len(data_2d))
    print("3d", len(data_3d))
    print("cfid3d", len(data_ml))
    print("ff", len(data_ff))
    assert (len(data_2d), len(data_3d), len(data_ml), len(data_ff)) == (
        1070,
        36099,
        35984,
        3291,
    )


# test_figshare_download()
