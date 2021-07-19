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


"""
if __name__ == "__main__":
    from jarvis.db.figshare import get_jid_data
    from jarvis.core.atoms import Atoms

    a = Atoms.from_dict(get_jid_data(jid="JVASP-667")["atoms"])
    cmd = "/home/knc6/Software/zeopp/zeo++-0.3/network"
    x, y, z = get_porosity(atoms=a, network_cmd=cmd)
    print(x, y, z)
"""
