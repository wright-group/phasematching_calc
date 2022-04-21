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
tkwat = 0.0001

thins = 100

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "water"
samp1.description = desc

samp1.load_layer(lay2file, tkwat, label="h2o")

# samp1.layers[0].suppress_absorbances()

# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [1800.0, 2700.0, 30000.0]
las.add_frequencies(arr1)
arr2 = [-5.0, 2.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")

var2 = np.linspace(2450.00, 2900.00, 91)[None, :]
var2a = np.linspace(2450.00, 2900.00, 91)
var1 = np.linspace(1300.0, 1900.0, 161)[:, None]
var1a = np.linspace(1300.0, 1900.0, 161)

ch1 = np.zeros([len(var1a), len(var2a)])
ch1p = np.zeros([len(var1a), len(var2a)])
ch1c = np.zeros([len(var1a), len(var2a)])
A1list = np.zeros([len(var1a), len(var2a)])
A2list = np.zeros([len(var1a), len(var2a)])
A3list = np.zeros([len(var1a), len(var2a)])

newtk = thins * tkwat

for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        Alist, Alistout = pc.phasematch.calculate_absorbances(samp1, las)

        ch1[m, n] = np.sqrt(Mlist[0])
        ch1p[m, n] = Mphase[0]
        A1list[m, n] = Alist[0][0]
        A2list[m, n] = Alist[0][1]
        A3list[m, n] = Alist[0][2]

for m in range(len(var1a)):
    for n in range(len(var2a)):
        Mconjsum = 0 + 0 * 1j
        for i in range(thins):
            Mphasedelta = (i + 1) * ch1p[m, n]
            Mconjdelta = np.cos(Mphasedelta) + np.sin(Mphasedelta) * 1j
            Mconj = (
                ch1[m, n]
                * Mconjdelta
                * 10 ** (-i * A1list[m, n])
                * 10 ** (-i * A2list[m, n])
                * 10 ** (-i * A3list[m, n])
                * tkwat
            )
            Mconjsum = Mconjsum + Mconj
        ch1c[m, n] = np.abs(Mconjsum) * np.abs(Mconjsum)
        pass

samp1.change_layer(1, thickness=newtk)

ch2 = np.zeros([len(var1a), len(var2a)])

for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(1, var1a[m])
        las.change_freq(2, var2a[n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        ch2[m, n] = Mlist[0] * newtk * newtk


A1list_t = np.zeros([len(var1a), len(var2a)])
A2list_t = np.zeros([len(var1a), len(var2a)])

A1list_t = A1list * thins
A2list_t = A2list * thins


data = wt.Data(name="example")
data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="DOVE", values=ch1)
data.create_channel(name="DOVESUM", values=ch1c)
data.create_channel(name="DOVETHICK", values=ch2)
data.create_channel(name="A1", values=(A1list_t))
data.create_channel(name="A2", values=(A2list_t))
data.transform("w1", "w2")
wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

wt.artists.quick2D(data, channel=2)
plt.show()

wt.artists.quick2D(data, channel=3)
plt.show()

wt.artists.quick2D(data, channel=4)
plt.show()

pass
