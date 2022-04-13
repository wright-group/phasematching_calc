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
samp1.load_layer(lay3file, tksapph, label="caf2fw")
samp1.load_layer(lay4file, tkacn, label="ACN")
samp1.load_layer(lay3file, tksapph, label="caf2bw")

# new Lasers object
las4=pc.Lasers.Lasers()
arr1=[3150.0,2200.0,17200.0]
las4.add_frequencies(arr1)
arr2=[6.0,-13.0,0.0]
las4.add_angles(arr2)
arr3=[1,-1,1]
las4.add_k_coeffs(arr3)
arr4=[1,1,1]
las4.add_pols(arr4)
las4.change_geometry("planar")

#out=pc.phasematch.m_plot(samp1,las4,2,2)

angl1=pc.phasematch.solve_angle(samp1,las4,2,2, isclose=True)
out=list(angl1)
print(out)

freq=pc.phasematch.solve_frequency(samp1,las4,2,3,20)
out=list(freq)
print(out[0])

las4.change_freq(3,out[0])

las4.change_freq(2,2190.0)
angle=pc.phasematch.solve_frequency(samp1,las4,2,3,20)
out2=list(angle)
print(out2[0])

las4.change_freq(3,out[0])
angle=pc.phasematch.solve_angle(samp1,las4,2,2)
out3=list(angle)
print(out3[0])



