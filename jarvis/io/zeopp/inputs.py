"""Module to run zeo++ package."""
import os
import tempfile


def get_porosity(atoms=None, network_cmd="./network", output_file=None):
    """Gete pore diameters using zeo++."""
    new_file, filename = tempfile.mkstemp()
    filename = filename + ".cif"
    atoms.write_cif(filename)
    if output_file is None:
        new_file, filename1 = tempfile.mkstemp()
        output_file = filename1
    cmd = network_cmd + " -ha -res " + output_file + " " + filename
    os.system(cmd)
    f = open(output_file, "r")
    lines = f.read().splitlines()
    f.close()
    largest_included_sphere = lines[0].split()[1]
    largest_free_sphere = lines[0].split()[2]
    largest_included_sphere_along_free_sphere_path = lines[0].split()[3]
    os.close(new_file)
    os.remove(filename)
    return (
        largest_included_sphere,
        largest_free_sphere,
        largest_included_sphere_along_free_sphere_path,
    )


def get_surface_area(
    atoms=None,
    network_cmd="/home/knc6/Software/zeopp/zeo++-0.3/network",
    probe_radius=1.2,
    chan_radius=None,
    num_samples=2000,
    output_file=None,
    prefix="ja_atoms",
):
    """Gete surface area using zeo++."""
    # See: http://www.zeoplusplus.org/examples.html
    if chan_radius is None:
        chan_radius = probe_radius
    filename = str(prefix) + ".cif"
    atoms.write_cif(filename)
    filename1 = str(prefix) + ".sa"
    if os.path.exists(filename1):
        os.remove(filename1)
    cmd = (
        network_cmd
        + " -ha -sa "
        + str(probe_radius)
        + " "
        + str(chan_radius)
        + " "
        + str(num_samples)
        + " "
        + filename
    )
    os.system(cmd)
    f = open(filename1, "r")
    lines = f.read().splitlines()
    f.close()
    ASA_m2cm3 = float(lines[0].split()[9])
    ASA_m2g = float(lines[0].split()[11])
    return ASA_m2cm3, ASA_m2g


"""
if __name__ == "__main__":
    from jarvis.db.figshare import get_jid_data
    from jarvis.core.atoms import Atoms

    a = Atoms.from_dict(get_jid_data(jid="JVASP-667")["atoms"])
    cmd = "/home/knc6/Software/zeopp/zeo++-0.3/network"
    x, y, z = get_porosity(atoms=a, network_cmd=cmd)
    print(x, y, z)
"""
