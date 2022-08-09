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

filepath = os.path.join(ROOT_DIR, "tests")

lay1file = os.path.join(filepath, "CaF2_Malitson.txt")
lay2file = os.path.join(filepath, "H2O_1.txt")
tksap = 0.02
tkwat = 0.0008

thins = 30
thick = thins * tkwat

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "water"
samp1.description = desc

samp1.load_layer(lay2file, tkwat, label="h2o")
# use lay1file for a limit where absorbance is zero

# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [1800.0, 2700.0, 30000.0]
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
A1list = np.zeros([len(var1a), len(var2a)])
A2list = np.zeros([len(var1a), len(var2a)])
A3list = np.zeros([len(var1a), len(var2a)])
w4 = np.zeros([len(var1a), len(var2a)])


# thin DOVE summation calculations...it needs to tabulate absorbance changes
# as each thin layer is added (A1list, A2list, A3list).   The phase change
# as a result in changing absorbances is found in ch1p.   The ch1 calculates
# a single thin layer m factor.

#  ---  Absorbance list in a single thin layer determined.
for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        w4t = (
            las.frequencies[0] * las.k_coeffs[0]
            + las.frequencies[1] * las.k_coeffs[1]
            + las.frequencies[2] * las.k_coeffs[2]
        )
        w4t2, a, n4t = samp1["layers"][0].estimate(w4t)
        Alist, Alistout = pc.phasematch.calculate_absorbances(samp1, las)
        A1list[m, n] = Alist[0][0]
        A2list[m, n] = Alist[0][1]
        A3list[m, n] = Alist[0][2]
        w4[m, n] = w4t

#  ---  following is a slow algorithm.
for m in range(len(var1a)):
    for n in range(len(var2a)):
        Mconjsum = 0.000 + 0.000 * 1j
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        w4t = w4[m, n]
        for i in range(thins):
            if i == 0:
                Mphaseprev = 0.000
            else:
                Mphaseprev = Mphase[0]
            E1power = np.sqrt(10 ** (-i * A1list[m, n] / 2.00))  # 2.00 converts I/Io to E/Eo
            E2power = np.sqrt(10 ** (-i * A2list[m, n] / 2.00))
            E3power = np.sqrt(10 ** (-i * A3list[m, n] / 2.00))
            tktemp = i * thins
            Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
            Mphase[0] = Mphase[0] + Mphaseprev
            Mlisttemp = np.sqrt(Mlist[0]) * E1power * E2power * E3power
            # Phase differential is calculated backwards from the final layer
            Mphasedelta = (
                np.cos(w4t * (2 * np.pi) * (thick - tktemp) + Mphase[0])
                + np.sin(w4t * (2 * np.pi) * (thick - tktemp) + Mphase[0]) * 1j
            )
            # if i ==(thins-1):
            #  Mphasedelta=0.000
            Mconjtemp = Mlisttemp * (Mphasedelta)
            Mconjsum = Mconjsum + Mconjtemp * tkwat
        ch1[m, n] = np.abs(Mconjsum) * np.abs(Mconjsum)


# thick Dove calculations
samp1.change_layer(1, thickness=thick)
ch2 = np.zeros([len(var1a), len(var2a)])

for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        ch2[m, n] = Mlist[0] * thick * thick

# A1list_t = np.zeros([len(var1a), len(var2a)])
# A2list_t = np.zeros([len(var1a), len(var2a)])
# A1list_t = A1list * thins
# A2list_t = A2list * thins

data = wt.Data(name="example")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="DOVESUM", values=ch1)
data.create_channel(name="DOVETHICK", values=ch2)
# data.create_channel(name="A1", values=(A1list_t))
# data.create_channel(name="A2", values=(A2list_t))
data.transform("w1", "w2")

wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

pass
