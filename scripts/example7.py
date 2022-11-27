import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *


# filepath = os.path.join(ROOT_DIR, "tests")
filepath = os.path.join(os.getcwd(), "tests")
lay1file = os.path.join(filepath, "sapphire1.txt")
lay2file = os.path.join(filepath, "H2O_1.txt")
tksap = 0.02
tkwat = 0.01


# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "FWM cell with fw sapphire, sample water, and bw sapphire"
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
ch1a = np.zeros([len(var1), len(var2a)])
ch2 = np.zeros([len(var1), len(var2a)])
ch2a = np.zeros([len(var1), len(var2a)])
ch3 = np.zeros([len(var1), len(var2a)])
ch3a = np.zeros([len(var1), len(var2a)])

for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        Alist, Alistout = pc.phasematch.calculate_absorbances(samp1, las)
        Mlista = pc.phasematch.apply_absorbances(Mlist, Alist, Alistout)
        Mlistb = pc.phasematch.apply_trans(Mlista, Tdict)
        samp1.change_layer(2, thickness=0.0001)
        Mlist1a, Mphase1a, tklist1a, Tdict1a = pc.phasematch.m_calc(samp1, las)
        ch1[m, n] = Mlist[1]
        ch1a[m, n] = Mlist1a[1]
        samp1.change_layer(2, thickness=tkwat)

vec2 = [1, 1, 1]
las.add_k_coeffs(vec2)

for m in range(len(var1)):
    for n in range(len(var2a)):
        las.change_freq(1, var1[m])
        las.change_freq(2, var2a[n])
        Mlist2, Mphase, tklist2, Tlist2 = pc.phasematch.m_calc(samp1, las)
        ch2[m, n] = Mlist2[1]
        samp1.change_layer(2, thickness=0.0001)
        Mlist2a, Mphase2a, tklist2a, Tdict2a = pc.phasematch.m_calc(samp1, las)
        ch2a[m, n] = Mlist2a[1]
        samp1.change_layer(2, thickness=tkwat)

ch3 = ch1 / ch2
ch3a = ch1a / ch2a


data = wt.Data(name="example")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="DOVE", values=ch1)
data.create_channel(name="TSF", values=ch2)
data.create_channel(name="DOVE_TSF_RATIO", values=ch3)
data.create_channel(name="DOVE_TSF_RATIO_thinfilm", values=ch3a)
data.transform("w1", "w2")
wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

wt.artists.quick2D(data, channel=2)
plt.show()

wt.artists.quick2D(data, channel=3)
plt.show()  # should be 1 for all data points or very close to it
pass
