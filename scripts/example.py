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


filepath = os.path.join(ROOT_DIR, "tests")

lay1file = os.path.join(filepath, "CaF2_Malitson.txt")
lay2file = os.path.join(filepath, "H2O_1.txt")

tkcaf2 = 0.02  # cm
tkwater = 0.01  # cm

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample(description="FWM Cell")

samp1.load_layer(lay1file, tkcaf2, label="caf2fw")
samp1.load_layer(lay2file, tkwater, label="water")
samp1.load_layer(lay1file, tkcaf2, label="caf2bw")

# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [1800.0, 2700.0, 30000.0]
las.add_frequencies(arr1)
arr2 = [20.0, 7.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

# Test scripts show how the above objects can be saved and loaded.

# single point Mcalc check
Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)


# angle estimation for laser 1 in layer 2 with frequency 1750 cm-1
# using default geometry
# angle=pc.phasematch.SolveAngle(samp1,las,2,1,1750)

# frequency estimate for a laser 1 in layer 2 to allow for phasematching
# with the angle in air shown in the Las object above.
# freq=pc.phasematch.SolveFrequency(samp1, las, 2, 1)


# A method for creating a 2D array of "Mcalcs" and converting into a
# WrightTools data object for use in various simulations.
var1 = np.linspace(2450.00, 2900.00, 46)[:, None]
var2 = np.linspace(1300.0, 1900.0, 61)[None, :]
var2a = np.linspace(1300.0, 1900.0, 61)

ch1 = np.zeros([len(var1), len(var2a)])
for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tlist = pc.phasematch.m_calc(samp1, las)
        ch1[m, n] = np.abs(Mlist[1])


data = wt.Data(name="FWM cell water w/CaF2 planar DOVE")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="Mfactor", values=ch1)
data.transform("w1", "w2")
wt.artists.quick2D(data)
plt.show()

pass
