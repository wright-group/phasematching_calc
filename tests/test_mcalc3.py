import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *


filepath=os.path.join(ROOT_DIR, 'tests')

lay3file=os.path.join(filepath, 'CaF2_Malitson.txt')
lay4file=os.path.join(filepath, 'CH3CN_paste_1.txt')
lay5file=os.path.join(filepath, 'Silicon_shk_ir.txt')

tkcaf2=0.02 #cm
tkacn=0.01 #cm
tksi=0.02 

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay5file, tksi, label="siliconfw")
samp1.loadlayer(lay4file, tkacn, label="ACN")
samp1.loadlayer(lay3file, tkcaf2, label="caf2bw")

# new Lasers object
las4=pc.Lasers.Lasers()
arr1=[3150.0,2200.0,5000.0]
las4.addfrequencies(arr1)
arr2=[20.0,15.0,0.0]
las4.addangles(arr2)
arr3=[1,-1,1]
las4.addkcoeffs(arr3)
arr4=[1,1,1]
las4.addpolarizations(arr4)
las4.changegeometry("planar")

angle=pc.phasematch.SolveAngle(samp1,las4,2,2,2200)
out=angle
print(out)

freq=pc.phasematch.SolveFrequency(samp1,las4,2,3,30)
out=freq
print(out)