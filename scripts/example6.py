import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *



filepath=os.path.join(ROOT_DIR, 'tests')


#new IsoSample: sapphire: ACN: sapphire
lay3file=os.path.join(filepath, 'CaF2_Malitson.txt')
lay4file=os.path.join(filepath, 'CH3CN_paste_1.txt')

tksapph=0.02 #cm
tkacn=0.01 #cm

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay3file, tksapph, label="caf2fw")
samp1.loadlayer(lay4file, tkacn, label="ACN")
samp1.loadlayer(lay3file, tksapph, label="caf2bw")

# new Lasers object
las4=pc.Lasers.Lasers()
arr1=[3150.0,2200.0,25000.0]
las4.addfrequencies(arr1)
arr2=[6.0,-15.0,0.0]
las4.addangles(arr2)
arr3=[1,-1,1]
las4.addkcoeffs(arr3)
arr4=[1,1,1]
las4.addpolarizations(arr4)
las4.changegeometry("planar")


freq=pc.phasematch.SolveFrequency(samp1,las4,2,3,20)
out=list(freq)
print(out[0])

las4.changefreq(3,out[0])

las4.changefreq(2,2190.0)
angle=pc.phasematch.SolveFrequency(samp1,las4,2,3,20)
out2=list(angle)
print(out2[0])

las4.changefreq(3,out[0])
angle=pc.phasematch.SolveAngle(samp1,las4,2,2)
out3=list(angle)
print(out3[0])



