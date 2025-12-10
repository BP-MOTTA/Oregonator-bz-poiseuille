# -*- coding: utf-8 -*-
"""
Created on Sun Feb 13 23:20:37 2022

@author: BP_motta
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from scipy.interpolate import griddata
archivos= glob.glob(f"*.dat")
posi=[[],[],[]]
tf=[]
for filename in archivos:
    t=int(filename[0:4])
    tf.append(t)
    nombre="t="+ str(t)[0:6]
    x=np.arange(1,800.1,0.10)
    y=np.arange(1,21.1,0.10)
    z=np.loadtxt(filename,usecols=[0,1,2],unpack=True)
    zi=griddata((z[0],z[1]),z[2],(x[None,:], y[:,None]))
    fig, axs = plt.subplots(2)
    zi=griddata((z[0],z[1]),z[2],(x[None,:], y[:,None]))
    axs[0].imshow(zi,origin="lower",vmin=0,aspect="auto",extent=[1,800,1,21])
    b=np.where(z[1]==5)
    xl=z[0,b].tolist()
    zl=z[2,b].tolist()
    axs[1].plot(xl[0],zl[0])
    fig.suptitle(nombre)
    noarchi=nombre+".jpg"
    plt.savefig(noarchi)
    