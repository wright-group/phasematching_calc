import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
""" Following example shows how the phasematching (M) factors can be calculated in a 2_dimensional
manner and converted to a WrightTools data object.  Other methods can be employed such as use
in a WrightSim object or in other simulations.  The other methods found in the phasematching script
are shown as well."""

#NOTE it is important to note that Mcalc does not multiply by thickness^2 currently.  This
# may be added later but may slightly limit the capability of the calculator as then l^4 dependencies
# such as cascades will not be estimatable.

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
arr1=[1800.0,2700.0,30000.0]
las.addfrequencies(arr1)
arr2=[10.0,-10.0, 0.0]
las.addangles(arr2)
arr3=[-1,1,1]
las.addkcoeffs(arr3)
arr4=[1,1,1]
las.addpolarizations(arr4)
las.changegeometry("planar")

#Test scripts show how the above objects can be saved and loaded.

# single point Mcalc check
Mlist,tklist,Tdict=pc.phasematch.Mcalc(samp1,las)

# angle estimation for laser 1 in layer 2 with frequency 1750 cm-1
# using default geometry
#angle=pc.phasematch.SolveAngle(samp1,las,2,1,1750)

# frequency estimate for a laser 1 in layer 2 to allow for phasematching
# with the angle in air shown in the Las object above.
#freq=pc.phasematch.SolveFrequency(samp1, las, 2, 1)


# A method for creating a 2D array of "Mcalcs" and converting into a
# WrightTools data object for use in various simulations.
var1=np.linspace(2500.00,2900.00,41)[:,None]
var2=np.linspace(1650.0,1900.0,36)[None, :]
var2a=np.linspace(1650.0,1900.0,36)

ch1= np.zeros([len(var1), len(var2a)])
for m in range(len(var1)):
    for n in range(len(var2a)):
        las.changefreq(1,var1[m])
        las.changefreq(2,var2a[n])
        Mlist,tklist,Tlist=pc.phasematch.Mcalc(samp1,las)
        ch1[m,n]=(np.abs(Mlist[2])**2)*tkwater**2  #final multiplier to get the M factor into a result 
                            # where values can be readily understood


data=wt.Data(name="example")
data.create_variable(name="w1", units="wn", values= var1)
data.create_variable(name="w2", units="wn", values= var2)
data.create_channel(name='Mfactor', values=ch1)
data.transform("w1","w2")
wt.artists.quick2D(data)
plt.show()

pass