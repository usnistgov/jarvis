import subprocess
from subprocess import Popen, PIPE
import sys
import os

dsfile2 = open("/home/smandal/wien_struc_init.py", "r")
content2 = []
nnodes = 1
nprocs = 10
jidlist = ["JVASP-1002"]

for lines in dsfile2:
    content2.append(lines)
cwd = "/oasis/scratch/comet/smandal/temp_project/Silicon"
for id in jidlist:
    dirname2 = str(id)  ##+str("_DMFT")
    f1L = str(os.getcwd()) + str("/") + str(id)
    if not os.path.exists(dirname2):
        os.makedirs(dirname2)

    os.chdir(dirname2)
    dir_file2 = open("dmft.py", "w")
    for lines in content2:
        dir_file2.write(lines)
    line = (
        str("prepare_wien_input(jid=") + str("'") + str(id) + str("'") + str(")") + "\n"
    )
    dir_file2.write(line)
    dir_file2.close()
    ref1 = open("REFERENCE", "w")
    line = str(id)
    ref1.write(line)
    ref1.close()
    ref2 = open("FUNCTIONAL", "w")
    line = str("PBE")
    ref2.write(line)
    ref2.close()
    f = open("submit_job", "w")
    line = str("#!/bin/bash") + "\n"
    f.write(line)
    line = str("#SBATCH -J ") + str(id) + "\n"
    f.write(line)
    line = str("#SBATCH -t 00:30:00") + "\n"
    f.write(line)
    line = str("#SBATCH --partition=debug") + "\n"
    f.write(line)
    line = str("#SBATCH -o test.log") + "\n"
    f.write(line)
    line = str("#SBATCH -A TG-DMR180109") + "\n"
    f.write(line)
    line = str("#SBATCH -N ") + str(nnodes) + "\n"
    # f.write(line)
    line = str("#SBATCH --ntasks-per-node=") + str(nprocs) + "\n"
    ##        f.write(line)
    ##        line=str("#SBATCH --mem=")+str(mem)+'\n'
    f.write(line)
    dir = str(os.getcwd())
    line = str("cd ") + str(cwd) + str("/") + str(id) + "\n"
    f.write(line)
    line = (
        str(
            "module load intel/2013_sp1.2.144 mkl/11.1.2.144 mvapich2_ib/2.1 fftw/3.3.8   "
        )
        + "\n"
    )
    ##line=str("module load intel/intel-2018.03 gsl/2.4  fftw/3.3.8-mpi   python/2.7.15 ")+'\n'
    f.write(line)
    line = (
        str(
            "export PATH=/home/smandal/miniconda2/envs/my_jarvis/bin:/home/smandal/miniconda2/bin:/home/smandal/miniconda2/condabin:/home/smandal/bini:/home/smandal/W2k/w2k_14:/home/smandal/gap2c:/home/smandal/bini:/home/smandal/W2k/w2k_14:/home/smandal/gap2c:/opt/gnu/gcc/bin:/opt/gnu/bin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/sbin:/opt/sdsc/bin:/opt/sdsc/sbin:/opt/ibutils/bin:/usr/java/latest/bin:/opt/pdsh/bin:/opt/rocks/bin:/opt/rocks/sbin:/home/smandal/W2k/w2k_14:/home/smandal/W2k/w2k_14/SRC_structeditor/bin:/home/smandal/W2k/w2k_14/SRC_IRelast/script-elastic:.:/home/smandal/W2k/w2k_14:.:/home/smandal/bin:/home/smandal/W2k/w2k_14:/home/smandal/W2k/w2k_14/SRC_structeditor/bin:/home/smandal/W2k/w2k_14/SRC_IRelast/script-elastic:.:/home/smandal/W2k/w2k_14:."
        )
        + "\n"
    )

    f.write(line)
    line = str("source activate my_jarvis") + "\n"
    f.write(line)
    line = str("export WIENROOT=/hom2/smandal/W2k/w2k_14/") + "\n"
    f.write(line)
    line = str("export PATH=$PATH:$WIENROOT:$WIEN_DMFT_ROOT") + "\n"
    f.write(line)
    ##line=str("export PYTHONPATH=$PYTHONPATH:.:/home1/sm2159/.local/lib/python2.7/site-packages/:$WIEN_DMFT_ROOT")+'\n'
    ##f.write(line)
    line = str("export MODULEPATH=$MODULEPATH:/home/smandal/modulefiles") + "\n"
    f.write(line)
    ##line=str(" export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/software/intel/ps-xe-2018.03/compilers_and_libraries_2018.3.222/linux/mkl/lib/intel64:/software/fftw/3.3.8-mpi/lib:/home1/sm2159/.local/lib/python2.7/site-packages/:/software/gsl/2.4/lib:$PATH")+'\n'
    ##f.write(line)
    ##line=str("export PATH=.:$PATH")+'\n'
    ##f.write(line)
    # line=str(" mpi_prefix='/software/intel/ps-xe-2018.03/compilers_and_libraries_2018.3.222/linux/mpi/intel64/bin/mpirun -n $SLURM_NTASKS'")+'\n'
    # f.write(line)
    line = str("export SLURM_NODEFILE=`generate_pbs_nodefile`") + "\n"
    f.write(line)
    line = str("cat $SLURM_NODEFILE > JOBFILE") + "\n"
    f.write(line)
    line = str("echo 'granularity:1' > ./.machines") + "\n"
    f.write(line)
    tmp = str('{print "1:"$1":3"}')
    line = str("awk '") + str(tmp) + str("' JOBFILE >> ./.machines") + "\n"
    f.write(line)
    line = str("echo 'extrafine:3'   >> ./.machines") + "\n"
    f.write(line)
    line = (
        str("/home/smandal/miniconda2/envs/my_jarvis/bin/python ")
        + str(cwd)
        + str("/")
        + str(id)
        + str("/dmft.py>out")
        + "\n"
    )
    f.write(line)
    f.close()
    with open("job.out", "w") as f:
        p = subprocess.Popen(
            ["sbatch", "submit_job"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = p.communicate()
        job_id = stdout  # .rstrip('\n').split()[-1]
        print("stdout,stderr", stdout, stderr)
        # job_id = str(stdout.split('Your job')[1].split(' ')[1])
        ##f.write(job_id)
    os.chdir("../")
