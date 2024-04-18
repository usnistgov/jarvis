from jarvis.core.composition import Composition


def test_comp():
    comp = {"Li": 2, "O": 4}
    cc = Composition(comp)

    assert (
        cc.prototype_new,
        cc.formula,
        cc.reduced_formula,
        round(cc.weight, 4),
        cc.to_dict(),
    ) == ("A2B", "Li2O4", "LiO2", 77.8796, comp)
    # ) == ("AB2", "Li2O4", "LiO2", 77.8796, comp)
    assert cc.prototype == "AB2"
    c = Composition.from_string("Al2O3Al5Co6O1")
    td = c.to_dict()
    fd = Composition.from_dict(td)
    assert c.formula == "Al7Co6O4"
    arr = c.atomic_fraction_array
