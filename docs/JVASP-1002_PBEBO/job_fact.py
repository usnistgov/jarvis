from jarvis.tasks.vasp.vasp import JobFactory
from jarvis.db.jsonutils import loadjson
d=loadjson("job_fact.json")
v=JobFactory.from_dict(d)
v.all_optb88vdw_calcs()
