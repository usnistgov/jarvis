"""Run multiple jobs."""
from jarvis.tasks.qe.super import SuperCond
from jarvis.core.utils import get_factors
from jarvis.core.atoms import Atoms
from jarvis.db.figshare import get_jid_data
from jarvis.core.kpoints import Kpoints3D
from jarvis.tasks.queue_jobs import Queue
from jarvis.db.jsonutils import dumpjson
import os
from jarvis.analysis.structure.spacegroup import Spacegroup3D
import glob


jids = ["JVASP-816"]
run_dir = "/home/knc6/Software/qe/jarvis/jarvis/tasks/qe"
run_dir = "/rk2/knc6/SuperCon/QE_RUN2"
jids1 = [
    "JVASP-88846",
    "JVASP-987",
    "JVASP-981",
    "JVASP-25388",
    "JVASP-984",
    "JVASP-834",
    "JVASP-14604",
    "JVASP-837",
    "JVASP-840",
    "JVASP-25379",
    "JVASP-25144",
    "JVASP-14744",
    "JVASP-890",
    "JVASP-888",
]
jids = jids1 + [
    "JVASP-4406",
    "JVASP-19821",
    "JVASP-14622",
    "JVASP-969",
    "JVASP-972",
    "JVASP-25254",
    "JVASP-25407",
    "JVASP-961",
    "JVASP-958",
    "JVASP-963",
    "JVASP-14832",
    "JVASP-966",
    "JVASP-25125",
    "JVASP-802",
    "JVASP-25273",
    "JVASP-25167",
    "JVASP-919",
    "JVASP-25114",
    "JVASP-922",
    "JVASP-949",
    "JVASP-95268",
    "JVASP-79561",
    "JVASP-1056",
    "JVASP-14612",
    "JVASP-102277",
    "JVASP-943",
    "JVASP-931",
    "JVASP-934",
    "JVASP-937",
    "JVASP-21193",
    "JVASP-946",
    "JVASP-25142",
    "JVASP-828",
    "JVASP-33718",
    "JVASP-1011",
    "JVASP-25250",
    "JVASP-25213",
    "JVASP-1002",
    "JVASP-14601",
    "JVASP-14812",
    "JVASP-14837",
    "JVASP-996",
    "JVASP-993",
    "JVASP-7804",
    "JVASP-858",
    "JVASP-25104",
    "JVASP-25180",
    "JVASP-852",
    "JVASP-25248",
    "JVASP-1035",
    "JVASP-864",
    "JVASP-861",
    "JVASP-867",
    "JVASP-910",
    "JVASP-25117",
    "JVASP-25337",
    "JVASP-916",
    "JVASP-1026",
    "JVASP-1029",
    "JVASP-25210",
    "JVASP-1017",
    "JVASP-1020",
    "JVASP-1014",
    "JVASP-21197",
    "JVASP-870",
    "JVASP-895",
    "JVASP-14725",
    "JVASP-1050",
    "JVASP-810",
    "JVASP-14606",
    "JVASP-901",
    "JVASP-816",
    "JVASP-14603",
    "JVASP-819",
    "JVASP-825",
    "JVASP-898",
    "JVASP-21195",
]

subd = (
    []
)  # ['JVASP-837', 'JVASP-25379', 'JVASP-972', 'JVASP-834', 'JVASP-958', 'JVASP-984', 'JVASP-25273', 'JVASP-961', 'JVASP-25125', 'JVASP-963', 'JVASP-25167', 'JVASP-966', 'JVASP-14622', 'JVASP-14744', 'JVASP-14832', 'JVASP-890', 'JVASP-802', 'JVASP-14604', 'JVASP-25407', 'JVASP-888', 'JVASP-987', 'JVASP-840', 'JVASP-25254', 'JVASP-969', 'JVASP-88846', 'JVASP-25144', 'JVASP-25388', 'JVASP-981']


# jids=['JVASP-834']
jids = ["JVASP-19821", "JVASP-14622"]
jids = ["JVASP-4406"]  # ,"JVASP-19821",'JVASP-14622']
jids = ["JVASP-19821"]
jids = ["JVASP-816", "JVASP-934", "JVASP-961", "JVASP-19821", "JVASP-4406"]
jids = ["JVASP-816"]
jids = ["JVASP-816", "JVASP-934", "JVASP-961", "JVASP-19821", "JVASP-4406"]
qe_cmd = "/cluster/deb9/bin/mpirun /cluster/bin/pw.x"
qe_cmd = "/cluster/bin/pw.x"
qe_cmd = "/home/knc6/Software/qe/q-e/bin/pw.x"

##########################################################################################
qe_cmd = "/users/knc6/anaconda2/envs/qe38/bin/pw.x"
qe_cmd = " /cluster/deb9/bin/mpirun -n 16 /cluster_sw/bin/pw.x"
jids = ["JVASP-19821", "JVASP-4406"]
run_dir = "/working/knc6/Sup"
jids = [
    "JVASP-25210",
    "JVASP-1029",
    "JVASP-1056",
    "JVASP-898",
    "JVASP-14832",
    "JVASP-993",
    "JVASP-895",
    "JVASP-25248",
    "JVASP-1050",
    "JVASP-25180",
    "JVASP-946",
    "JVASP-837",
    "JVASP-21193",
    "JVASP-25250",
    "JVASP-14603",
    "JVASP-1020",
    "JVASP-870",
    "JVASP-25379",
    "JVASP-943",
    "JVASP-828",
    "JVASP-864",
    "JVASP-25167",
    "JVASP-1014",
    "JVASP-14744",
    "JVASP-888",
    "JVASP-958",
    "JVASP-25104",
    "JVASP-14837",
    "JVASP-969",
    "JVASP-931",
    "JVASP-919",
    "JVASP-949",
    "JVASP-996",
    "JVASP-972",
    "JVASP-1017",
    "JVASP-916",
    "JVASP-890",
    "JVASP-1026",
    "JVASP-14612",
    "JVASP-819",
    "JVASP-25273",
    "JVASP-95268",
    "JVASP-910",
    "JVASP-963",
    "JVASP-966",
    "JVASP-25117",
    "JVASP-802",
    "JVASP-25407",
    "JVASP-102277",
    "JVASP-1011",
    "JVASP-852",
    "JVASP-981",
    "JVASP-825",
    "JVASP-14606",
    "JVASP-984",
    "JVASP-14604",
    "JVASP-21195",
    "JVASP-810",
    "JVASP-937",
    "JVASP-901",
    "JVASP-7804",
    "JVASP-25254",
    "JVASP-858",
    "JVASP-79561",
    "JVASP-861",
    "JVASP-25142",
    "JVASP-834",
    "JVASP-987",
    "JVASP-922",
    "JVASP-25114",
    "JVASP-840",
    "JVASP-25125",
    "JVASP-25213",
    "JVASP-25337",
    "JVASP-33718",
    "JVASP-14622",
    "JVASP-14601",
    "JVASP-25388",
    "JVASP-14725",
    "JVASP-25144",
    "JVASP-867",
]


jids = [
    "JVASP-19821",
    "JVASP-19684",
    "JVASP-19668",
    "JVASP-11981",
    "JVASP-14960",
    "JVASP-19767",
    "JVASP-19679",
    "JVASP-934",
    "JVASP-20620",
    "JVASP-961",
    "JVASP-816",
]
jids = ["JVASP-15938", "JVASP-81987"]
jids = [
    "JVASP-1020",
    "JVASP-14741",
    "JVASP-14837",
    "JVASP-14601",
    "JVASP-1014",
    "JVASP-898",
    "JVASP-14492",
    "JVASP-25337",
]
jids = ["JVASP-4406"]
# screened lambda
jids = [
    "JVASP-14961",
    "JVASP-20082",
    "JVASP-18513",
    "JVASP-243",
    "JVASP-19933",
    "JVASP-14984",
    "JVASP-63394",
    "JVASP-1378",
    "JVASP-12066",
    "JVASP-18711",
    "JVASP-19",
    "JVASP-14802",
    "JVASP-14680",
    "JVASP-15801",
    "JVASP-28357",
    "JVASP-359",
    "JVASP-93855",
    "JVASP-14926",
]
# a15_mgb2
jids = [
    "JVASP-14961",
    "JVASP-4406",
    "JVASP-63045",
    "JVASP-63844",
    "JVASP-63770",
    "JVASP-14853",
    "JVASP-85341",
    "JVASP-18660",
    "JVASP-15815",
    "JVASP-14795",
    "JVASP-15816",
    "JVASP-19926",
    "JVASP-15080",
    "JVASP-95772",
    "JVASP-14798",
    "JVASP-92307",
    "JVASP-1378",
    "JVASP-93926",
    "JVASP-93930",
    "JVASP-20096",
    "JVASP-1735",
    "JVASP-20071",
    "JVASP-94297",
    "JVASP-14940",
    "JVASP-14647",
    "JVASP-19907",
    "JVASP-14867",
    "JVASP-15015",
    "JVASP-19933",
    "JVASP-100603",
    "JVASP-20082",
    "JVASP-103658",
    "JVASP-19934",
    "JVASP-100962",
    "JVASP-4861",
    "JVASP-19836",
    "JVASP-19821",
    "JVASP-14680",
    "JVASP-18598",
    "JVASP-14726",
    "JVASP-14954",
    "JVASP-19827",
    "JVASP-18931",
    "JVASP-13993",
    "JVASP-14905",
    "JVASP-18932",
    "JVASP-18640",
    "JVASP-14736",
    "JVASP-19921",
    "JVASP-7807",
    "JVASP-8723",
    "JVASP-18652",
    "JVASP-18654",
    "JVASP-14407",
    "JVASP-4771",
    "JVASP-7734",
    "JVASP-14097",
    "JVASP-7745",
    "JVASP-7942",
    "JVASP-56594",
    "JVASP-1151",
    "JVASP-78264",
    "JVASP-78528",
    "JVASP-78530",
    "JVASP-78278",
    "JVASP-78492",
    "JVASP-78582",
    "JVASP-78493",
    "JVASP-78934",
    "JVASP-78500",
    "JVASP-78810",
    "JVASP-78472",
    "JVASP-78813",
    "JVASP-78814",
    "JVASP-78818",
    "JVASP-78480",
    "JVASP-36220",
    "JVASP-78512",
    "JVASP-36221",
    "JVASP-78662",
    "JVASP-36264",
    "JVASP-78887",
    "JVASP-78888",
    "JVASP-36232",
    "JVASP-36375",
    "JVASP-36410",
    "JVASP-36379",
    "JVASP-58273",
    "JVASP-56195",
    "JVASP-36302",
    "JVASP-56933",
    "JVASP-36434",
    "JVASP-36206",
    "JVASP-35084",
    "JVASP-20237",
    "JVASP-14254",
    "JVASP-20433",
    "JVASP-20743",
    "JVASP-20374",
    "JVASP-19596",
    "JVASP-16481",
    "JVASP-39386",
    "JVASP-15916",
    "JVASP-19628",
    "JVASP-16532",
    "JVASP-19643",
    "JVASP-19774",
    "JVASP-20272",
    "JVASP-19735",
    "JVASP-37611",
    "JVASP-16711",
    "JVASP-13755",
    "JVASP-35211",
    "JVASP-4325",
    "JVASP-37996",
    "JVASP-20645",
    "JVASP-14499",
    "JVASP-20624",
    "JVASP-17781",
    "JVASP-14544",
    "JVASP-4319",
    "JVASP-20565",
    "JVASP-20596",
    "JVASP-19959",
    "JVASP-16466",
    "JVASP-16470",
    "JVASP-16472",
    "JVASP-16475",
    "JVASP-19805",
    "JVASP-16504",
    "JVASP-16505",
    "JVASP-16546",
    "JVASP-19961",
    "JVASP-16548",
    "JVASP-16550",
    "JVASP-16552",
    "JVASP-16556",
    "JVASP-19970",
    "JVASP-19671",
    "JVASP-17234",
    "JVASP-19685",
    "JVASP-54597",
    "JVASP-15901",
    "JVASP-19730",
    "JVASP-19641",
    "JVASP-19745",
    "JVASP-19723",
    "JVASP-20506",
    "JVASP-20538",
    "JVASP-4678",
    "JVASP-15897",
    "JVASP-4373",
    "JVASP-17526",
    "JVASP-16869",
    "JVASP-36236",
    "JVASP-36253",
    "JVASP-54925",
    "JVASP-114151",
    "JVASP-115780",
    "JVASP-114249",
    "JVASP-116158",
    "JVASP-114455",
    "JVASP-116202",
    "JVASP-116311",
    "JVASP-114601",
    "JVASP-115135",
    "JVASP-115181",
    "JVASP-115191",
    "JVASP-115482",
    "JVASP-113425",
    "JVASP-113491",
    "JVASP-113497",
    "JVASP-113503",
    "JVASP-113511",
    "JVASP-113526",
    "JVASP-113663",
    "JVASP-114023",
    "JVASP-114063",
    "JVASP-114089",
    "JVASP-118884",
    "JVASP-118886",
    "JVASP-118898",
    "JVASP-118132",
    "JVASP-120216",
    "JVASP-118601",
    "JVASP-120242",
    "JVASP-117640",
    "JVASP-117679",
    "JVASP-117704",
    "JVASP-120316",
    "JVASP-120336",
    "JVASP-117863",
    "JVASP-120382",
    "JVASP-117949",
    "JVASP-120955",
    "JVASP-121186",
    "JVASP-121206",
    "JVASP-121285",
    "JVASP-91617",
    "JVASP-14960",
    "JVASP-20086",
    "JVASP-90822",
    "JVASP-91708",
    "JVASP-91577",
    "JVASP-15777",
    "JVASP-90875",
    "JVASP-14703",
    "JVASP-15003",
    "JVASP-92231",
    "JVASP-19860",
    "JVASP-91560",
    "JVASP-91607",
    "JVASP-18952",
    "JVASP-15008",
    "JVASP-14642",
    "JVASP-18419",
    "JVASP-18790",
    "JVASP-18793",
    "JVASP-18418",
    "JVASP-15046",
    "JVASP-19834",
    "JVASP-14722",
    "JVASP-19837",
    "JVASP-11981",
    "JVASP-19917",
    "JVASP-18771",
    "JVASP-7743",
    "JVASP-56595",
    "JVASP-56193",
    "JVASP-20275",
    "JVASP-20371",
    "JVASP-20246",
    "JVASP-20434",
    "JVASP-56974",
    "JVASP-19630",
    "JVASP-19632",
    "JVASP-16454",
    "JVASP-19703",
    "JVASP-19777",
    "JVASP-18246",
    "JVASP-14771",
    "JVASP-14513",
    "JVASP-14769",
    "JVASP-18149",
    "JVASP-17794",
    "JVASP-18048",
    "JVASP-14506",
    "JVASP-20626",
    "JVASP-18146",
    "JVASP-20562",
    "JVASP-20516",
    "JVASP-19800",
    "JVASP-19668",
    "JVASP-19976",
    "JVASP-19684",
    "JVASP-20544",
    "JVASP-19719",
    "JVASP-19785",
    "JVASP-19755",
    "JVASP-19650",
    "JVASP-15938",
    "JVASP-19765",
    "JVASP-19749",
    "JVASP-19769",
    "JVASP-16715",
    "JVASP-110513",
]


def non_prime_kpoints(kpts=[]):
    mem = []
    for i in kpts:
        facts = get_factors(i)
        if len(facts) == 1:
            val = i + 1
        else:
            val = i
        mem.append(val)
    return mem


# jids=['JVASP-19821', 'JVASP-14960', 'JVASP-19679', 'JVASP-934', 'JVASP-961', 'JVASP-1014', 'JVASP-14601', 'JVASP-816', 'JVASP-20620', 'JVASP-19684', 'JVASP-19767', 'JVASP-14741', 'JVASP-14837', 'JVASP-898', 'JVASP-1020', 'JVASP-11981', 'JVASP-14492', 'JVASP-15938', 'JVASP-19668', 'JVASP-25337', 'JVASP-81987']


def write_qejob(pyname="job.py", job_json=""):
    """Write template job.py with VaspJob.to_dict() job.json."""
    f = open(pyname, "w")
    f.write("from jarvis.tasks.qe.super import SuperCond\n")
    f.write("from jarvis.db.jsonutils import loadjson\n")
    f.write('d=loadjson("' + str(job_json) + '")\n')
    f.write("v=SuperCond.from_dict(d)\n")
    f.write("v.runjob()\n")
    f.close()


submit_job = True
submitted = []
for i in glob.glob("JVASP*"):
    submitted.append(i.split("_")[0])
print(len(submitted), submitted)

import sys

# sys.exit()
for i in jids:
    if i not in submitted:
        try:
            print("jid", i)
            dir_name = os.path.join(run_dir, i + "_SUPER")
            if not os.path.exists(dir_name):
                os.makedirs(dir_name)
            os.chdir(dir_name)
            dat = get_jid_data(jid=i, dataset="dft_3d")
            a_atoms = Atoms.from_dict(dat["atoms"])
            atoms = Spacegroup3D(a_atoms).refined_atoms.get_primitive_atoms
            print(atoms)
            kp = Kpoints3D().automatic_length_mesh(
                # lattice_mat=atoms.lattice_mat,
                # length=10
                lattice_mat=atoms.lattice_mat,
                length=dat["kpoint_length_unit"],
            )
            kpts = kp._kpoints[0]
            kpts = non_prime_kpoints(kpts)
            kp = Kpoints3D(kpoints=[kpts])
            print("kpts", kpts)

            nq1 = get_factors(kpts[0])[0]
            nq2 = get_factors(kpts[1])[0]
            nq3 = get_factors(kpts[2])[0]
            qp = Kpoints3D(kpoints=[[nq1, nq2, nq3]])

            sup = SuperCond(atoms=atoms, kp=kp, qp=qp, qe_cmd=qe_cmd).to_dict()
            dumpjson(data=sup, filename="sup.json")
            write_qejob(job_json=os.path.abspath("sup.json"))
            path = (
                "echo hello"
                + " \nsource ~/anaconda2/envs/my_jarvis/bin/activate my_jarvis \npython "
                + os.getcwd()
                + "/job.py"
            )
            if submit_job:
                directory = os.getcwd()
                tmp_name = str(i) + "_submit_job"
                submit_cmd = ["sbatch", tmp_name]
                # For Slurm clusters
                Queue.slurm(
                    job_line=path,
                    filename=tmp_name,
                    jobname=i,
                    queue="rack2,rack3,rack4,rack4e,rack5,rack6,rack6i",
                    directory=directory,
                    submit_cmd=submit_cmd,
                    walltime="100:00:00",
                )

                """
     Queue.pbs(
        job_line=path,
        jobname=i,
        directory=os.getcwd(),
        submit_cmd=["qsub", "-q","highmem","submit_job"],
     )
     """

            os.chdir("..")
        except:
            pass
