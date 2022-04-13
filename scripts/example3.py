import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *



filepath=os.path.join(ROOT_DIR, 'tests')

lay1file=os.path.join(filepath, 'CH3CN_paste_1.txt')
lay2file=os.path.join(filepath, 'sapphire1.txt')


tksap=0.02
tkacn=0.01 

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.load_layer(lay1file, tksap, label="sapphire")
samp1.load_layer(lay2file, tkacn, label="acn")
samp1.load_layer(lay1file, tksap, label="sapphire")


#generation of a Lasers object.
las=pc.Lasers.Lasers()
arr1=[1800.0,2700.0,18400.0]
las.add_frequencies(arr1)
arr2=[8.0,-7.0, 0.0]
las.add_angles(arr2)
arr3=[-1,1,1]
las.add_k_coeffs(arr3)
arr4=[1,1,1]
las.add_pols(arr4)
las.change_geometry("planar")

var1=np.linspace(2600.00,3200.00,61)[:,None]
var2=np.linspace(1600.0,2200.0,61)[None, :]
var2a=np.linspace(1600.0,2200.0,61)

ch1= np.zeros([len(var1), len(var2a)])
mold=int(0)
for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1,var1[m])
        las.change_freq(2,var2a[n])
        if (m==0 & n==0):
            angleair2=pc.phasematch.solve_angle(samp1,las,2,1,isclose=False)
            angletemp=list(angleair2)[1]
            ch1[m,n]=list(angleair2)[1]
            las.change_angle(2,list(angleair2)[1])  
        elif (mold==m):
            angleair2=pc.phasematch.solve_angle(samp1,las,2,1,isclose=True)
            ch1[m,n]=(list(angleair2)[0]) 
            las.change_angle(2,list(angleair2)[0])           
        else:
            las.change_angle(2,angletemp) 
            angleair2=pc.phasematch.solve_angle(samp1,las,2,1,isclose=True)
            ch1[m,n]=list(angleair2)[0]
            mold=m
            angletemp=list(angleair2)[0]
            las.change_angle(2,list(angleair2)[0])  


data=wt.Data(name="angle check for w2")
data.create_variable(name="w1", units="wn", values= var1)
data.create_variable(name="w2", units="wn", values= var2)
data.create_channel(name='angleforw2', values=ch1)
data.transform("w2","w1")
wt.artists.quick2D(data)
plt.show()

pass