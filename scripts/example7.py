import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *



filepath=os.path.join(ROOT_DIR, 'tests')

lay1file=os.path.join(filepath, 'sapphire1.txt')
lay2file=os.path.join(filepath, "H2O_1.txt")
tksap=0.02 
tkwat=0.01


# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="sapphwatersapph"
samp1.description=desc
samp1.loadlayer(lay1file, tksap, label="saphfw")
samp1.loadlayer(lay2file, tkwat, label="h2o")
samp1.loadlayer(lay1file, tksap, label="saphfw")


#generation of a Lasers object.
las=pc.Lasers.Lasers()
arr1=[1800.0,2700.0,30000.0]
las.addfrequencies(arr1)
arr2=[-18.0,8.0, 0.0]
las.addangles(arr2)
arr3=[-1,1,1]
las.addkcoeffs(arr3)
arr4=[1,1,1]
las.addpolarizations(arr4)
las.changegeometry("planar")


#angle1=pc.phasematch.SolveAngle(samp1,las,1,1,frequency=1800.0)
#print(list(angle1))
#las.changeangle(1,list(angle1)[0])

var1=np.linspace(2450.00,2900.00,91)[:,None]
var2=np.linspace(1300.0,1900.0,161)[None, :]
var2a=np.linspace(1300.0,1900.0,161)

ch1= np.zeros([len(var1), len(var2a)])
ch2=np.zeros([len(var1), len(var2a)])
ch3=np.zeros([len(var1), len(var2a)])

chartin,chartout=pc.phasematch.calculatedeltats(samp1, las)

for m in range(len(var1)):
    for n in range(len(var2a)):
        las.changefreq(1,var1[m])
        las.changefreq(2,var2a[n])
        Mlist,tklist,Tdict=pc.phasematch.Mcalc(samp1,las)
        Alist, Alistout=pc.phasematch.calculateabsorbances(samp1,las) 
        Mlist1a=pc.phasematch.applyabsorbances(Mlist,Alist,Alistout)
        Mlist1b=pc.phasematch.applyfresneltrans(Mlist1a, Tdict)
        ch1[m,n]=Mlist[1]  

vec2=[1,1,1]
las.addkcoeffs(vec2)
pass

for m in range(len(var1)):
    for n in range(len(var2a)):
        las.changefreq(1,var1[m])
        las.changefreq(2,var2a[n])
        Mlist2,tklist2,Tlist2=pc.phasematch.Mcalc(samp1,las)
        ch2[m,n]=Mlist2[1]  

ch3=ch1/ch2

data=wt.Data(name="example")
data.create_variable(name="w1", units="wn", values= var1)
data.create_variable(name="w2", units="wn", values= var2)
data.create_channel(name='DOVE', values=ch1)
data.create_channel(name='TSF', values=ch2)
data.create_channel(name="DOVE_TSF_RATIO", values=ch3)
data.transform("w1","w2")
wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

wt.artists.quick2D(data, channel=2)
plt.show()

pass