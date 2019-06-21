import numpy as np
#from lib import view, traj_view
import matplotlib.pyplot as plt
from IPython.display import Image
import json as json 
from ipywidgets import interact, fixed
from elements import ELEMENTS
"""
This is a collection of functions needed to manage the interactive
plots and the parsing of the databases in the ML exercise
  AA + TS 
"""

def parser(d):
    mats=[]
    Ks=[]
    for m in d:
        pos=m['poscar']
        lines=pos.splitlines()
        a=np.array([float(x) for x in lines[2].split()])
        b=np.array([float(x) for x in  lines[3].split()])
        c=np.array([float(x) for x in lines[4].split()])
        #if c[2]< 0: c=-c
        species=lines[5].split()
        spe_count=[int(x) for x in lines[6].split()]
        i=8
        p=[]
        mat={'species':[], 'positions':[],'cell':[a,b,c]}
        for j in range(sum(spe_count)):
            mat['positions'].append(sum([float(x)*np.array(mat['cell'][k]) for k,x in enumerate(lines[i+j].split()[:3])]))
            mat['species'].append(ELEMENTS[lines[i+j].split()[3]].number)
        mats.append(mat)
        #
        Ks.append(m["K_VRH"])
    return mats, np.array(Ks)
