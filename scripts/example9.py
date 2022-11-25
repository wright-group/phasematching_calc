import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *

"""The following example slices a single, thick layer into multiple very small layers and compares
the sum (with phases) to the single layer to see if they are nearly equal.  Errors in how the phase
or M factor are calculated may manifest as difference and so this example is actually a test for how well
the calculation works.  Calculation does not perform absorbance of the output vector as the conditions
(water, visible output) expected little absorbance there.  This example tests the validity of the phase
portion of the m_calc method."""

# filepath = os.path.join(ROOT_DIR, "tests")
filepath = os.path.join(os.getcwd(), "tests")
# filepath=os.getcwd()
lay1file = os.path.join(filepath, "sapphire1.txt")
lay2file = os.path.join(filepath, "H2O_1.txt")
tksap = 0.02
tkwat = 0.01

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "sapphire fw, water sample, sapphire bw"
samp1.description = desc
samp1.load_layer(lay1file, tksap, label="saphfw")
samp1.load_layer(lay2file, tkwat, label="h2o")
samp1.load_layer(lay1file, tksap, label="saphfw")


# generation of a Lasers object.
las = pc.Lasers.Lasers()
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

data = wt.Data(name="example")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="angleoutx(deg)", values=ch1)

data.transform("w1", "w2")
wt.artists.quick2D(data, channel=0)
plt.show()


pass
