import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *
import time

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
arr2 = [-13.0, 10.0, 0.0]
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
test1 = np.zeros([len(var1a), len(var2a)])

mold = int(0)
for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[n])
        las.change_freq(2, var2a[m])
        if (m == 0) & (n == 0):
            angleair2, amount = (pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=False))
            angletemp = list(angleair2)[0]  # this needs to solve for remainder to work
            if np.any(list(angleair2)):
                ch1[m, n] = angletemp
                las.change_angle(1, angletemp)
        elif mold == m:
            angleair2a, amt = (pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount))
            if np.any(list(angleair2)):
                ch1[m, n] = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])
        else:
            las.change_angle(1, angletemp)
            angleair2, amt = (pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount))
            mold = m
            if np.any(list(angleair2)):
                ch1[m, n] = list(angleair2)[0]
                angletemp = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])

data = wt.Data(name="angle for lower frequency input, opp side")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="angleforw1", values=ch1)
data.transform("w2", "w1")
data.channels[0].signed=True

wt.artists.quick2D(data)
plt.show()


for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[n])
        las.change_freq(2, var2a[m])
        if (m == 0) & (n == 0):
            angleair2, amount = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=False)
            angletemp = list(angleair2)[1]  # this needs to solve for remainder to work
            if np.any(list(angleair2)):
                ch1[m, n] = angletemp
                las.change_angle(1, angletemp)
        elif mold == m:
            angleair2, amt = (pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount))
            if np.any(list(angleair2)):
                ch1[m, n] = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])
        else:
            las.change_angle(1, angletemp)
            angleair2, amt = (pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True, amt=amount))
            mold = m
            if np.any(list(angleair2)):
                ch1[m, n] = list(angleair2)[0]
                angletemp = list(angleair2)[0]
                las.change_angle(1, list(angleair2)[0])


for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[n])
        las.change_freq(2, var2a[m])
        las.change_angle(1, ch1[m, n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        test1[m, n] = Mlist[1]


data2 = wt.Data(name="angle for lower frequency beam, same side")
data2.create_variable(name="w1", units="wn", values=var1)
data2.create_variable(name="w2", units="wn", values=var2)
data2.create_channel(name="angleforw1", values=ch1)
data2.channels[0].signed=True
data2.create_channel(name="test", values=test1)
data2.transform("w2", "w1")
wt.artists.quick2D(data2, channel=0)
plt.show()

wt.artists.quick2D(data2, channel=1)
plt.show()
