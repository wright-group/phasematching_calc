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

#Other test scripts show how the above objects can be saved and loaded.
# single point Mcalc check
Mlist,tklist,Tdict=pc.phasematch.Mcalc(samp1,las)
testoutput=np.abs(Mlist[1])**2*tkwater**2
print(testoutput)


# angle estimation for laser 1 in layer 2 with frequency 2600 cm-1
# using boxcars..test for all reals
angle=pc.phasematch.SolveAngle(samp1,las,2,1,2600)
out=angle
print(out)


# frequency estimate for a laser 1 in layer 2 to allow for phasematching
# with the angle in air shown in the Las object above...test for all reals >0
freq=pc.phasematch.SolveFrequency(samp1, las, 2, 1)
out=freq
print(out)


#new lasers object:  unphasematchable geometry TSF
las2=pc.Lasers.Lasers()
arr1=[2700.0,1800.0,25000.0]
las2.addfrequencies(arr1)
arr2=[15.0,5.0,10.0]
las2.addangles(arr2)
arr3=[1,1,1]
las2.addkcoeffs(arr3)
arr4=[1,1,1]
las.addpolarizations(arr4)
las.changegeometry()

angle=pc.phasematch.SolveAngle(samp1,las2,2,1,2600)
out=angle
print(out)


#new lasers object:  planar geometry DOVE
las3=pc.Lasers.Lasers()
arr1=[2700.0,1500.0,38000.0]
las3.addfrequencies(arr1)
arr2=[15.0,5.0,10.0]
las3.addangles(arr2)
arr3=[1,-1,1]
las3.addkcoeffs(arr3)
arr4=[1,1,1]
las3.addpolarizations(arr4)
las3.changegeometry("planar")

angle=pc.phasematch.SolveAngle(samp1,las3,2,2,1600)
out=angle
print(out)

freq=pc.phasematch.SolveFrequency(samp1,las2,2,3)
out=freq
print(out)


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

freq=pc.phasematch.SolveFrequency(samp1,las4,2,3)
out=freq
print(out)

lay3file=os.path.join(filepath, 'sapphire1.txt')
lay4file=os.path.join(filepath, 'CH3CN_paste_1.txt')
lay5file=os.path.join(filepath, 'Silicon_shk_ir.txt')

tksapph=0.02 #cm
tkacn=0.01 #cm
tksi=0.02 

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay5file, tksi, label="siliconfw")
samp1.loadlayer(lay4file, tkacn, label="ACN")
samp1.loadlayer(lay3file, tksapph, label="sapphirebw")

# new Lasers object
las4=pc.Lasers.Lasers()
arr1=[3150.0,2200.0,5000.0]
las4.addfrequencies(arr1)
arr2=[80.0,15.0,0.0]
las4.addangles(arr2)
arr3=[1,-1,1]
las4.addkcoeffs(arr3)
arr4=[1,1,1]
las4.addpolarizations(arr4)
las4.changegeometry("planar")

angle=pc.phasematch.SolveAngle(samp1,las4,2,2,2200)
out=angle
print(out)

freq=pc.phasematch.SolveFrequency(samp1,las4,2,3)
out=freq
print(out)