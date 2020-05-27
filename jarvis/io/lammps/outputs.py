"""
Functions to analze LAMMPS output
"""

import sys


def analyze_log(log="log.lammps"):

    """ 
    Analyzes log.lammps file,
    Please note, the output format heavily depends on the input file
    A generic inpu is taken here
    
    Args:
    
        log: path to log.lammps file
        
    Returns:
    
          en: energy/atom
          
          press: pressure
          
          toten: total energy
          
          cij: elastic constants
    """
    en = 0
    press = 0
    c11 = 0
    c22 = 0
    c33 = 0
    c44 = 0
    c55 = 0
    c66 = 0
    c12 = 0
    c13 = 0
    c23 = 0
    c14 = 0
    c15 = 0
    c16 = 0
    c14 = 0
    c24 = 0
    c25 = 0
    c26 = 0
    c34 = 0
    c35 = 0
    c36 = 0
    c45 = 0
    c46 = 0
    c56 = 0
    try:
        logfile = open(log, "r")
        lines = logfile.read().splitlines()
        for i, line in enumerate(lines):
            if "Loop time of" in line:
                toten = float(lines[i - 1].split()[12])
                press = float(lines[i - 1].split()[2])
                press = float(press) * 0.0001
                en = float(lines[i - 1].split()[12]) / float(lines[i - 1].split()[17])
                break
        logfile.close()
    except:
        pass
    try:
        logfile = open(log, "r")
        lines = logfile.read().splitlines()
        for i, line in enumerate(lines):
            if 'print "Elastic Constant C11all = ${C11all} ${cunits}"' in line:
                c11 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C22all = ${C22all} ${cunits}"' in line:
                c22 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C33all = ${C33all} ${cunits}"' in line:
                c33 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C12all = ${C12all} ${cunits}"' in line:
                c12 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C13all = ${C13all} ${cunits}"' in line:
                c13 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C23all = ${C23all} ${cunits}"' in line:
                c23 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C44all = ${C44all} ${cunits}"' in line:
                c44 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C55all = ${C55all} ${cunits}"' in line:
                c55 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C66all = ${C66all} ${cunits}"' in line:
                c66 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C14all = ${C14all} ${cunits}"' in line:
                c14 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C16all = ${C16all} ${cunits}"' in line:
                c16 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C24all = ${C24all} ${cunits}"' in line:
                c24 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C25all = ${C25all} ${cunits}"' in line:
                c25 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C26all = ${C26all} ${cunits}"' in line:
                c26 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C34all = ${C34all} ${cunits}"' in line:
                c34 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C35all = ${C35all} ${cunits}"' in line:
                c35 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C36all = ${C36all} ${cunits}"' in line:
                c36 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C45all = ${C45all} ${cunits}"' in line:
                c45 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C46all = ${C46all} ${cunits}"' in line:
                c46 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
            if 'print "Elastic Constant C56all = ${C56all} ${cunits}"' in line:
                c56 = ((str((lines[i + 1])).split("=")[1]).split("GPa"))[0]
        logfile.close()
    except:
        pass
    return (
        round(en, 2),
        round(press, 2),
        float(toten),
        round(float(c11), 1),
        round(float(c22), 1),
        round(float(c33), 1),
        round(float(c12), 1),
        round(float(c13), 1),
        round(float(c23), 1),
        round(float(c44), 1),
        round(float(c55), 1),
        round(float(c66), 1),
        round(float(c14), 1),
        round(float(c16), 1),
        round(float(c24), 1),
        round(float(c25), 1),
        round(float(c26), 1),
        round(float(c34), 1),
        round(float(c35), 1),
        round(float(c36), 1),
        round(float(c45), 1),
        round(float(c46), 1),
        round(float(c56), 1),
    )

"""
if __name__ == "__main__":
    lg = "/rk2/knc6/JARVIS-FF/ALLOY8/Mishin_updated-Ni-Al-Co-2013.eam.alloy_nist/bulk@mp-134_fold/bulk@mp-134/log.lammps"
    x = analyze_log(lg)
    print(x)
"""
