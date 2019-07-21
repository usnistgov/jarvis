"""
Integration of JARVIS-FF/DFT/ML .
"""
import os
from pathlib import Path

try:
    # dir = os.path.dirname(os.path.realpath('__file__'))
    home = str(Path.home())
    #env = str(home) + str("/env_variables")
    env = str("env_variables")
    # env=str('env_variables'))
    f = open(env, "r")
    lines = f.read().splitlines()
    f.close()
    for i in lines:
        tmp = i.split("=")
        if len(tmp) > 1:
            os.environ[tmp[0]] = tmp[1]
except:
    pass
