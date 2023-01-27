import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from sympy import *

"""The following example calculates the changes in launch angle of DOVE FWM signal for an x planar geometry input with
fixed input angles but varying frequencies for w1 and w2.  This can be useful for making certain the full output is
captured by collection optics for any change of inputs, and can extend to other changes such as angles and w3."""

filepath = os.path.join(os.getcwd(), "tests")
lay1file = os.path.join(filepath, "sapphire1.txt")
lay2file = os.path.join(filepath, "H2O_1.txt")

tksap = 0.02
tkwat = 0.01

# generation of a IsoSample
# samp1 = pc.IsoSample.IsoSample()
samp1 = pc.IsoSample()
desc = "sapphire fw, water sample, sapphire bw"
samp1.description = desc
samp1.load_layer(lay1file, tksap, label="saphfw")
samp1.load_layer(lay2file, tkwat, label="h2o")
samp1.load_layer(lay1file, tksap, label="saphfw")


# generation of a Lasers object.
# las = pc.Lasers.Lasers()
las = pc.Lasers()
arr1 = [1800.0, 2700.0, 16000.0]
las.add_frequencies(arr1)
arr2 = [5.0, -2.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

var2 = np.linspace(2150.00, 3650.00, 151)[None, :]
var2a = np.linspace(2150.00, 3650.00, 151)
var1 = np.linspace(1200.0, 1900.0, 71)[:, None]
var1a = np.linspace(1200.0, 1900.0, 71)

ch1 = np.zeros([len(var1a), len(var2a)])


for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        angleoutx, angleouty = pc.phasematch.launchangle(samp1, las)
        ch1[m, n] = angleoutx

data = wt.Data(name="example 9")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="angleoutx(deg)", values=ch1)

data.transform("w1", "w2")
wt.artists.quick2D(data, channel=0)
plt.show()


pass
