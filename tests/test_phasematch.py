import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import S, FiniteSet, Interval, oo



filepath=os.path.join(ROOT_DIR, 'tests')

lay1file=os.path.join(filepath, 'CaF2_Malitson.txt')
lay2file=os.path.join(filepath, 'H2O_1.txt')

tkcaf2=0.02 #cm
tkwater=0.01 #cm

# generation of a IsoSample
samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay1file, tkcaf2, label="caf2fw")
samp1.loadlayer(lay2file, tkwater, label="water")
samp1.loadlayer(lay1file, tkcaf2, label="caf2bw")

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
testoutput=np.abs(Mlist[1])
assert (np.isclose(testoutput, 0.02978456))

# angle estimation for laser 1 in layer 2 with frequency 2600 cm-1
# using boxcars..test for all angles
angle=pc.phasematch.SolveAngle(samp1,las,2,1,2600)
assert (np.isclose(float(angle.end), float(77.9637423)))   # angle is an Interval, Interval.end is a Sympy function

# frequency estimate for a laser 1 in layer 2 to allow for phasematching
# with the angle in air shown in the Las object above...test for all reals >0
freq=pc.phasematch.SolveFrequency(samp1, las, 2, 1, 10)
out=freq
assert out==Interval(0,oo)   #this may be updated if critical angles will be calculated at high freq in future


Alist_in, Alist_out=pc.phasematch.calculateabsorbances(samp1,las)
assert (np.isclose(Alist_in[1][0], 1.180156))
assert (np.isclose(Alist_in[1][1], 2.973818))

tlist_in, tlist_out=pc.phasematch.calculatedeltats(samp1,las)
assert (np.isclose(tlist_in[1][0], 1413.241232))
assert (np.isclose(tlist_out[1], 1406.194641))

Mlistout=pc.phasematch.applyabsorbances(Mlist,Alist_in,Alist_out)
print(Mlistout)

#new lasers object:  unphasematchable geometry TSF
las2=pc.Lasers.Lasers()
arr1=[2700.0,1800.0,25000.0]
las2.addfrequencies(arr1)
arr2=[15.0,5.0,10.0]
las2.addangles(arr2)
arr3=[1,1,1]
las2.addkcoeffs(arr3)
arr4=[1,1,1]
las2.addpolarizations(arr4)
las2.changegeometry("boxcars")

angle=pc.phasematch.SolveAngle(samp1,las2,2,1,2600)
out=angle
assert out==S.EmptySet
