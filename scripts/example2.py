import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *


filepath = os.path.join(ROOT_DIR, "tests")

lay1file = os.path.join(filepath, "CaF2_Malitson.txt")

tkcaf2 = 0.03

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "caf2window300um"
samp1.description = desc
samp1.load_layer(lay1file, tkcaf2, label="caf2")


# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [1800.0, 2700.0, 18400.0]
las.add_frequencies(arr1)
arr2 = [8.0, 8.0, 8.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("boxcars")


var1 = np.linspace(2600.00, 3200.00, 61)[:, None]
var2 = np.linspace(1600.0, 2200.0, 61)[None, :]
var2a = np.linspace(1600.0, 2200.0, 61)

ch1 = np.zeros([len(var1), len(var2a)])
for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tlist = pc.phasematch.m_calc(samp1, las)
        ch1[m, n] = np.abs(Mlist[0])


data = wt.Data(name="CaF2 300 micron boxcars DOVE")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="Mfactor", values=ch1)
data.transform("w2", "w1")
wt.artists.quick2D(data)
plt.show()

pass
