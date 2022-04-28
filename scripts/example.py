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

"""Example1.  Generation of 2D phasematching plot."""

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
arr2 = [22.3, -8.6, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

# A method for creating a 2D array of "Mcalcs" and converting into a
# WrightTools data object for use in various simulations.
var1 = np.linspace(2450.00, 2850.00, 41)[:, None]
var2 = np.linspace(1450.0, 2000.0, 56)[None, :]
var2a = np.linspace(1450.0, 2000.0, 56)

ch1 = np.zeros([len(var1), len(var2a)])
for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tlist = pc.phasematch.m_calc(samp1, las)
        ch1[m, n] = np.abs(Mlist[1])

data = wt.Data(name="FWM cell water caf2 planar DOVE")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="Mfactor", values=ch1)
data.transform("w1", "w2")
wt.artists.quick2D(data)
plt.show()
