from jarvis.core.composition import Composition


def test_comp():
    comp = {"Li": 2, "O": 4}
    cc = Composition(comp)

    assert (
        cc.prototype,
        cc.formula,
        cc.reduced_formula,
        round(cc.weight, 4),
        cc.to_dict(),
    ) == ("AB2", "Li2O4", "LiO2", 77.8796, comp)
    c = Composition.from_string("Al2O3Al5Co6O1")
    td = c.to_dict()
    fd = Composition.from_dict(td)
    assert c.formula == "Al7Co6O4"
<<<<<<< HEAD
    arr = c.atomic_fraction_array
=======
>>>>>>> f5d6302fe6b017f8080f53cbf4a6d0110b18cfb1
