import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *



filepath=os.path.join(ROOT_DIR, 'tests')

lay1file=os.path.join(filepath, 'sapphire1.txt')
lay2file=os.path.join(filepath, 'H2O_1.txt')

tksapph=0.02 #cm
tkwater=0.01 #cm

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay1file, tksapph, label="sapphirefw")
samp1.loadlayer(lay2file, tkwater, label="water")
samp1.loadlayer(lay1file, tksapph, label="sapphirebw")

#generation of a Lasers object.
las=pc.Lasers.Lasers()
arr1=[2700.0,1800.0,1800.0]
las.addfrequencies(arr1)
arr2=[10.0,10.0,10.0]
las.addangles(arr2)
arr3=[1,-1,1]
las.addkcoeffs(arr3)
arr4=[1,1,1]
las.addpolarizations(arr4)
las.changegeometry("boxcars")

#Test scripts show how the above objects can be saved and loaded.
# single point Mcalc check
Mlist,Tdict=pc.phasematch.Mcalc(samp1,las)
testoutput=np.abs(Mlist[1])**2*tkwater**2
print(testoutput)


# angle estimation for laser 1 in layer 2 with frequency 2600 cm-1
# using boxcars..test for all reals
angle=pc.phasematch.EstimateAngle(samp1,las,2,1,2600)
if (angle.compare(S.Reals) ==0):
   print ("angles: all real numbers")
   pass
else: 
    out=list(angle)
    print(out)


# frequency estimate for a laser 1 in layer 2 to allow for phasematching
# with the angle in air shown in the Las object above...test for all reals >0
freq=pc.phasematch.EstimateFrequency(samp1, las, 2, 1)
if (freq.compare(S.Reals) ==0):
   print ("freqs: all real numbers >0")
   pass
else: 
    out=list(freq)
    print(out)


