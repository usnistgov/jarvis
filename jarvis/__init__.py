"""
Integration of JARVIS-FF/DFT/ML with genetic algorithm.
"""
import os

# env=str(os.path.join(os.path.dirname(__file__),'env_variables'))
try:
    dir = os.path.join(os.path.dirname(__file__))
    env = os.path.join(dir, "env_variables")
    # env=str('env_variables'))
    f = open(env, "r")
    lines = f.read().splitlines()
    f.close()
    for i in lines:
        tmp = i.split("=")
        if len(tmp) > 1:
            # print (tmp[0],tmp[1])
            os.environ[tmp[0]] = tmp[1]
except:
    print("Cannot find env_variables in site-packages/path")
    pass
