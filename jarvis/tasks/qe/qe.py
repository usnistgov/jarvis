"""Module to run QE jobs."""

from jarvis.io.qe.inputs import QEinfile
from jarvis.io.qe.outputs import QEout
import os
from jarvis.db.jsonutils import loadjson, dumpjson
import subprocess
from jarvis.core.atoms import Atoms
from jarvis.core.kpoints import Kpoints3D


class QEjob(object):
    """Module to run generic QE job."""

    def __init__(
        self,
        atoms=None,
        kpoints=None,
        input_params={},
        qe_cmd="pw.x",
        jobname="test",
        psp_dir=None,
        url=None,
        output_file="qe.out",
        input_file="qe.in",
        stderr_file="std.err",
        psp_temp_name=None,
    ):
        """Intitialize class."""
        self.atoms = atoms
        self.kpoints = kpoints
        self.input_params = input_params
        self.qe_cmd = qe_cmd
        self.jobname = jobname
        self.psp_dir = psp_dir
        self.url = url
        self.psp_temp_name = psp_temp_name
        self.input_file = input_file
        self.output_file = output_file
        self.stderr_file = stderr_file
        self.qeinput = QEinfile(
            atoms=self.atoms,
            kpoints=self.kpoints,
            psp_dir=self.psp_dir,
            input_params=self.input_params,
            url=self.url,
            psp_temp_name=self.psp_temp_name,
        )

    def write_input(self):
        """Write inputs."""
        self.qeinput.write_file(self.input_file)

    def to_dict(self):
        """Get dictionary."""
        info = {}
        info["atoms"] = self.atoms.to_dict()

        info["kpoints"] = self.kpoints.to_dict()
        info["qe_cmd"] = self.qe_cmd
        info["psp_dir"] = self.psp_dir
        info["url"] = self.url
        return info

    @classmethod
    def from_dict(self, info={}):
        """Load from a dictionary."""
        return QEjob(
            atoms=Atoms.from_dict(info["atoms"]),
            kpoints=Kpoints3D.from_dict(info["kpoints"]),
            qe_cmd=info["qe_cmd"],
            psp_dir=info["psp_dir"],
            url=info["url"],
        )

    def runjob(self):
        """Run job and make or return a metadata file."""
        fname = self.jobname + ".json"
        info = {}
        if not os.path.exists(fname):
            self.write_input()
            self.run()
            info = self.parse_outputs()
            if info["job_done"]:
                dumpjson(data=info, filename=fname)
        else:
            info = loadjson(fname)
        return info

    def run(self):
        """Use subprocess to tun a job."""
        with open(self.output_file, "w") as f_std, open(
            self.stderr_file, "w", buffering=1
        ) as f_err:
            # use line buffering for stderr
            cmd = self.qe_cmd + "<" + self.input_file
            print("cmd", cmd)
            p = subprocess.Popen(cmd, shell=True, stdout=f_std, stderr=f_err)
            p.wait()
        return p

    def parse_outputs(self):
        """Parse outputs from .out or .xml files."""
        info = {}
        info["out_path"] = "na"
        info["xml_path_"] = "na"
        info["total_energy"] = "na"
        info["job_done"] = "na"
        out_path = os.path.abspath(self.output_file)
        if os.path.exists(out_path):
            try:
                info["out_path"] = out_path
                qe_out = QEout(filename=out_path)
                job_done = qe_out.job_done
                info["job_done"] = job_done
                total_energy = qe_out.get_total_energy()
                info["total_energy"] = total_energy
            except Exception as exp:
                print("Exception", exp)
                pass
        if (
            "control" in self.input_params
            and "prefix" in self.input_params["control"]
        ):
            print("HERE1")
            xml_path = os.path.abspath(
                os.path.join(
                    self.input_params["control"]["prefix"]
                    .strip('"')
                    .strip("'")
                    + ".save",
                    "data-file-schema.xml",
                )
            )
            print("HERE2", xml_path)
            if os.path.exists(xml_path):
                info["xml_path"] = xml_path
        return info
