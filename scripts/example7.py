import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *


filepath = os.path.join(ROOT_DIR, "tests")

lay1file = os.path.join(filepath, "sapphire1.txt")
lay2file = os.path.join(filepath, "H2O_1.txt")
tksap = 0.02
tkwat = 0.01


# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "sapphwatersapph"
samp1.description = desc
samp1.load_layer(lay1file, tksap, label="saphfw")
samp1.load_layer(lay2file, tkwat, label="h2o")
samp1.load_layer(lay1file, tksap, label="saphfw")


# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [1800.0, 2700.0, 30000.0]
las.add_frequencies(arr1)
arr2 = [-18.0, 8.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

var1 = np.linspace(2450.00, 2900.00, 91)[:, None]
var2 = np.linspace(1300.0, 1900.0, 161)[None, :]
var2a = np.linspace(1300.0, 1900.0, 161)

ch1 = np.zeros([len(var1), len(var2a)])
ch2 = np.zeros([len(var1), len(var2a)])
ch3 = np.zeros([len(var1), len(var2a)])

for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        Alist, Alistout = pc.phasematch.calculate_absorbances(samp1, las)
        Mlist1a = pc.phasematch.apply_absorbances(Mlist, Alist, Alistout)
        Mlist1b = pc.phasematch.apply_trans(Mlist1a, Tdict)
        ch1[m, n] = Mlist[1]

vec2 = [1, 1, 1]
las.add_k_coeffs(vec2)

for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist2, tklist2, Tlist2 = pc.phasematch.m_calc(samp1, las)
        ch2[m, n] = Mlist2[1]

ch3 = ch1 / ch2

data = wt.Data(name="example")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="DOVE", values=ch1)
data.create_channel(name="TSF", values=ch2)
data.create_channel(name="DOVE_TSF_RATIO", values=ch3)
data.transform("w1", "w2")
wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

wt.artists.quick2D(data, channel=2)
plt.show()

pass
