"""Modules for job submission."""

import os
import subprocess
from collections import OrderedDict


class Queue(object):
    """Construct HPC Job class."""

    def __init__(
        self,
        q_type="head_node",
        q_parameters={},
        job_sub_cmd=None,
        job_check_cmd="",
        job_id="",
    ):
        """Initialize the class."""
        self.q_type = q_type
        self.q_parameters = q_parameters
        self.job_sub_cmd = job_sub_cmd
        self.job_check_cmd = job_check_cmd
        self.job_id = job_id

    def head_node(self, submit_cmd=None):
        """Select if run on the head node, not recommended."""
        os.system(submit_cmd)

    def to_dict(self):
        """Convert class to a dictionary."""
        d = OrderedDict()
        d["q_type"] = self.q_type
        d["q_parameters"] = self.q_parameters
        d["job_sub_cmd"] = self.job_sub_cmd
        d["job_check_cmd"] = self.job_check_cmd
        d["job_id"] = self.job_id
        return d

    @classmethod
    def from_dict(self, d={}):
        """Load from a dictionary."""
        return Queue(
            q_type=d["q_type"],
            q_parameters=d["q_parameters"],
            job_sub_cmd=d["job_sub_cmd"],
            job_check_cmd=d["job_check_cmd"],
            job_id=d["job_id"],
        )

    @classmethod
    def pbs(
        self,
        filename="submit_job",
        shell="#!/bin/bash",
        nnodes=1,
        cores=16,
        walltime=None,
        queue=None,
        account=None,
        group_name=None,
        jobname="myJob",
        jobout="job.out",
        joberr="job.err",
        memory=None,
        email=None,
        pre_job_lines=None,
        directory=None,
        env=None,
        job_line="echo I am here",
        post_job_lines=None,
        submit_cmd=None,
    ):
        """Select if run using PBS script."""
        f = open(filename, "w")
        f.write("%s\n" % shell)
        f.write("#PBS -l nodes=%d:ppn=%d\n" % (nnodes, cores))
        f.write("#PBS -N %s\n" % jobname)
        f.write("#PBS -o %s\n" % jobout)
        f.write("#PBS -e %s\n" % joberr)

        if walltime is not None:
            if isinstance(walltime, str):
                f.write("#PBS -l walltime=%s\n" % walltime)
            else:
                ValueError("Provide walltime in a string format", walltime)
        if queue is not None:
            f.write("#PBS -q %s\n" % queue)
        if account is not None:
            f.write("#PBS -A %s\n" % account)
        if group_name is not None:
            f.write("#PBS -G %s\n" % group_name)
        if env is not None:
            f.write("#PBS -V %s\n" % env)
        if email is not None:
            f.write("#PBS -M %s\n" % email)
        if memory is not None:
            f.write("#PBS -l %s\n" % memory)
        if pre_job_lines is not None:
            f.write("%s\n" % pre_job_lines)
        if directory is not None:
            f.write("cd %s\n" % directory)
        if job_line is not None:
            f.write("%s\n" % job_line)
        if post_job_lines is not None:
            f.write("%s\n" % post_job_lines)
        f.close()

        if submit_cmd is not None:

            with open("jobid", "w") as f:
                p = subprocess.Popen(
                    submit_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                stdout, stderr = p.communicate()
                print("stdout,stderr", stdout, stderr)
                # job_id = str(stdout.split('Your job')[1].split(' ')[1])
                f.write(str(stdout))

        print("x")

    @classmethod
    def slurm(
        self,
        filename="submit_job",
        shell="#!/bin/bash",
        nnodes=1,
        cores=16,
        walltime=None,
        queue=None,
        account=None,
        group_name=None,
        jobname="myJob",
        jobout="job.out",
        joberr="job.err",
        memory=None,
        email=None,
        pre_job_lines=None,
        directory=None,
        env=None,
        job_line="echo I am here",
        post_job_lines=None,
        submit_cmd=None,
    ):
        """Select if run using SLURM script."""
        f = open(filename, "w")
        f.write("%s\n" % shell)
        f.write("#SBATCH --nodes=%d\n" % (nnodes))
        f.write("#SBATCH --ntasks-per-node=%d\n" % (cores))

        if walltime is not None:
            if isinstance(walltime, str):
                f.write("#SBATCH --time=%s\n" % walltime)
            else:
                ValueError("Provide walltime in a string format", walltime)

        if queue is not None:
            f.write("#SBATCH --partition=%s\n" % (queue))

        if account is not None:
            f.write("#SBATCH --account=%s\n" % (account))

        if memory is not None:
            f.write("#SBATCH --mem=%s\n" % (memory))

        f.write("#SBATCH --error=%s\n" % joberr)
        f.write("#SBATCH --output=%s\n" % jobout)
        if pre_job_lines is not None:
            f.write("%s\n" % pre_job_lines)
        if directory is not None:
            f.write("cd %s\n" % directory)
        if job_line is not None:
            f.write("%s\n" % job_line)
        if post_job_lines is not None:
            f.write("%s\n" % post_job_lines)
        f.close()

        if submit_cmd is not None:

            with open("jobid", "w") as f:
                p = subprocess.Popen(
                    submit_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                stdout, stderr = p.communicate()
                print("stdout,stderr", stdout, stderr)
                # job_id = str(stdout.split('Your job')[1].split(' ')[1])
                f.write(str(stdout))


"""
if __name__ == "__main__":
    Queue().pbs()
    Queue().slurm(filename="ll")
    Queue().head_node("echo ABC")
"""
