import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *
import time

"""Example 3.  Simulation of a map of expected angles to achieve phasematching in the liquid layer
of a multilayer cell.  Plots separate both possible solutions (negative and positive)."""

filepath = os.path.join(ROOT_DIR, "tests")

lay1file = os.path.join(filepath, "CH3CN_paste_1.txt")
lay2file = os.path.join(filepath, "sapphire1.txt")
lay3file = os.path.join(filepath, "CaF2_Malitson.txt")

tksap = 0.02
tkacn = 0.01
tkcaf2 = 0.02

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "FWM cell"
samp1.description = desc
# samp1.load_layer(lay1file, tksap, label="sapphire")
samp1.load_layer(lay3file, tkcaf2, label="caf2")
samp1.load_layer(lay1file, tkacn, label="acn")
samp1.load_layer(lay3file, tkcaf2, label="caf2bw")
# samp1.load_layer(lay1file, tksap, label="sapphire")


# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [2200.0, 3150.0, 17200.0]
las.add_frequencies(arr1)
arr2 = [18.0, 8.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

var1 = np.linspace(1600.0, 2200.0, 61)[None, :]
var1a = np.linspace(1600.0, 2200.0, 61)
var2 = np.linspace(2600.00, 3200.00, 61)[:, None]
var2a = np.linspace(1600.0, 2200.0, 61)

ch1 = np.zeros([len(var1a), len(var2a)])
ch2 = np.zeros([len(var1a), len(var2a)])
test1 = np.zeros([len(var1a), len(var2a)])
test2 = np.zeros([len(var1a), len(var2a)])

mold = int(0)
for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[n])
        las.change_freq(2, var2a[m])
        if (m == 0) & (n == 0):
            """The first data point calculates the angle in a slow method."""
            angleair2, amount = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=False)
            angletemp = list(angleair2)[0]  # this needs to solve for remainder to work
            if np.any(list(angleair2)):
                ch1[m, n] = angletemp
                las.change_angle(1, angletemp)
        elif mold == m:
            """Afterwards it proceeds with a solve using the faster method. Unfortunately, this
            method may skip to the other solution if conditions are (un)favorable.   (Un)favorable conditions
            include heavy oscillations and bad initial guess for the isclose value."""
            angleair2a, amt = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount)
            if np.any(list(angleair2)):
                ch1[m, n] = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])
        else:
            """This final step is testing whether it is better to use the original solve upon a
            new scanline or to stick with the recent solve.  Currently in place is to roll back to
            the original solve.  It then updates angletemp for the next scanline."""
            las.change_angle(1, angletemp)
            angleair2, amt = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount)
            mold = m
            if np.any(list(angleair2)):
                ch1[m, n] = list(angleair2)[0]
                angletemp = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])

data = wt.Data(name="angle solves")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)

time1 = time.time()
# Other solution.
for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[n])
        las.change_freq(2, var2a[m])
        if (m == 0) & (n == 0):
            angleair2, amount = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=False)
            angletemp = list(angleair2)[1]  # this needs to solve for remainder to work
            if np.any(list(angleair2)):
                ch2[m, n] = angletemp
                las.change_angle(1, angletemp)
        elif mold == m:
            angleair2, amt = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount)
            if np.any(list(angleair2)):
                ch2[m, n] = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])
        else:
            las.change_angle(1, angletemp)
            angleair2, amt = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount)
            mold = m
            if np.any(list(angleair2)):
                ch2[m, n] = list(angleair2)[0]
                angletemp = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])
time2 = time.time()

print(time2 - time1)

for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[n])
        las.change_freq(2, var2a[m])
        las.change_angle(1, ch1[m, n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        las.change_angle(1, ch2[m, n])
        Mlist2, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        test1[m, n] = -np.log10(Mlist[1])
        test2[m, n] = -np.log10(Mlist2[1])


data.create_channel(name="angleforw1_negative", values=ch1)
data.channels[0].signed = True
data.create_channel(name="angleforw1_positive", values=ch2)
data.channels[1].signed = True
data.create_channel(
    name="test1", values=test1
)  # Tests to see if all M factors calculated are good
data.create_channel(
    name="test2", values=test2
)  # Tests to see if all M factors calculated are good
data.transform("w2", "w1")
wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

wt.artists.quick2D(data, channel=2)
plt.show()

wt.artists.quick2D(data, channel=3)
plt.show()
