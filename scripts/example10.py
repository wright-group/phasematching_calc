import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *
import copy

"""Simulation of a map of expected angles to achieve phasematching in planar geometry in the liquid layer
of a multilayer sapphire cell.   Positive solutions for w2 only, then M factor plotted for a single angle.   Test of
the create_layer w/ mole fraction method.   Note that the mole fractioned values of n and alpha do not necssarily
sum correctly because mixtures can interact with each other to modify results."""


# filepath = os.path.join(ROOT_DIR, "tests")
filepath = os.path.join(os.getcwd(), "tests")
lay1afile = os.path.join(filepath, "CCl4_paste.txt")
lay2file = os.path.join(filepath, "MeOH_kozma_ApplSpec_paste.txt")
lay3file = os.path.join(filepath, "sapphire1.txt")

tkccl4 = 0.0127  # cm
tksap = 0.016  # cm

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "FWM cell with fw sapphire, sample ccl4:meoh, and bw sapphire"
samp1.description = desc
samp1.load_layer(lay3file, tksap, label="sapfw")

layfilelist = list()
molfraclist = [0.9, 0.10]
layfilelist.append(lay1afile)
layfilelist.append(lay2file)

samp1.create_layer(layfilelist, molfraclist, thickness=tkccl4, label="ccl4_meoh")
samp1.load_layer(lay3file, tksap, label="sapbw")

# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [2200.0, 3150.0, 18500.0]
las.add_frequencies(arr1)
arr2 = [13.0, -8.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

w1start = 1800.00
w1end = 2000.00
w2start = 2800.000
w2end = 3000.00

var1 = np.linspace(w1start, w1end, 41)[None, :]
var1a = np.linspace(w1start, w1end, 41)
var2 = np.linspace(w2start, w2end, 41)[:, None]
var2a = np.linspace(w2start, w2end, 41)

ch1 = np.zeros([len(var1a), len(var2a)])
ch2 = np.zeros([len(var1a), len(var2a)])
ch3 = np.zeros([len(var1a), len(var2a)])
test1 = np.zeros([len(var1a), len(var2a)])
test2 = np.zeros([len(var1a), len(var2a)])

mold = int(0)

data = wt.Data(name="angle solve positive")
data.create_variable(name="w1", units="wn", values=var1.T)
data.create_variable(name="w2", units="wn", values=var2.T)

#  positive solution.
for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
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
                ch2[m, n] = float("nan")
        else:
            las.change_angle(1, angletemp)
            angleair2, amt = pc.phasematch.solve_angle(samp1, las, 2, 1, isclose=True)
            mold = m
            if np.any(list(angleair2)):
                ch2[m, n] = list(angleair2)[-1]
                angletemp = list(angleair2)[-1]
                las.change_angle(1, list(angleair2)[-1])
            else:
                ch2[m, n] = float("nan")

# fixed angle M factor solver
angle1 = 13

las.change_angle(1, angle1)
las.change_angle(2, arr2[1])

for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tlist = pc.phasematch.m_calc(samp1, las)
        ch3[m, n] = np.abs(Mlist[1])


Mlist, Mphase, tklist, Tlist = pc.phasematch.m_calc(samp1, las)

data.transform("w2", "w1")

data.create_channel(name="angleforw1_positive", values=ch2)
data.channels[0].null = np.min(ch2)

data.create_channel(name=f"Mfactor at angle1={angle1} deg", values=ch3)
data.channels[1].null = 0
# data.channels[0].null = 0

wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

pass
