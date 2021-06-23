from jarvis.core.specie import (
    Specie,
    get_feats_hot_encoded,
    get_digitized_feats_hot_encoded,
    get_specie_data,
    atomic_numbers_to_symbols
)


def test_sp():
    el = Specie("Al")
    assert (
        el.Z,
        round(el.atomic_mass, 2),
        el.symbol,
        round(el.get_chgdescrp_arr[1], 2),
        round(el.get_descrp_arr[1], 2),
    ) == (13, 26.98, "Al", 12.17, 2792.11)
    dat = get_feats_hot_encoded()
    fat = get_digitized_feats_hot_encoded()
    x, y, z = get_specie_data()
    x=atomic_numbers_to_symbols()
