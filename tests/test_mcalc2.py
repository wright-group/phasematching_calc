import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *



filepath=os.path.join(ROOT_DIR, 'tests')


#new IsoSample: sapphire: ACN: sapphire
lay3file=os.path.join(filepath, 'sapphire1.txt')
lay4file=os.path.join(filepath, 'CH3CN_paste_1.txt')

tksapph=0.02 #cm
tkacn=0.01 #cm

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay3file, tksapph, label="sapphirefw")
samp1.loadlayer(lay4file, tkacn, label="ACN")
samp1.loadlayer(lay3file, tksapph, label="sapphirebw")

# new Lasers object
las4=pc.Lasers.Lasers()
arr1=[3150.0,2200.0,12500.0]
las4.addfrequencies(arr1)
arr2=[15.0,5.0,10.0]
las4.addangles(arr2)
arr3=[1,-1,1]
las4.addkcoeffs(arr3)
arr4=[1,1,1]
las4.addpolarizations(arr4)
las4.changegeometry("planar")

angle=pc.phasematch.SolveAngle(samp1,las4,2,2,2200)
out=angle
print(out)

freq=pc.phasematch.SolveFrequency(samp1,las4,2,3,20)
out=freq
print(out)
