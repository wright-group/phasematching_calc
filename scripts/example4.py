import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *


filepath = os.path.join(ROOT_DIR, "tests")

lay1file = os.path.join(filepath, "CH3CN_paste_1.txt")
lay2file = os.path.join(filepath, "sapphire1.txt")


tksap = 0.02
tkacn = 0.01

# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "FWM cell"
samp1.description = desc
samp1.load_layer(lay1file, tksap, label="sapphire")
samp1.load_layer(lay2file, tkacn, label="acn")
samp1.load_layer(lay1file, tksap, label="sapphire")


# generation of a Lasers object.
las = pc.Lasers.Lasers()
arr1 = [1800.0, 2700.0, 24000.0]
las.add_frequencies(arr1)
arr2 = [12.0, -7.0, 0.0]
las.add_angles(arr2)
arr3 = [-1, 1, 1]
las.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las.add_pols(arr4)
las.change_geometry("planar")


var1 = np.linspace(2600.00, 3200.00, 61)[:, None]
var1a = np.linspace(2600.00, 3200.00, 61)
var2 = np.linspace(1600.0, 2200.0, 61)[None, :]
var2a = np.linspace(1600.0, 2200.0, 61)

ch1 = np.zeros([len(var1a), len(var2a)])
test1 = np.zeros([len(var1a), len(var2a)])

mold = int(0)
for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(2, var1a[m])
        las.change_freq(1, var2a[n])
        if (m == 0) & (n == 0):
            """The first data point calculates the frequency in a slower method."""
            freqsolve, amount = pc.phasematch.solve_frequency(samp1, las, 2, 3, isclose=False)
            freqtemp = list(freqsolve)[0]  # this needs to solve for remainder to work
            if np.any(list(freqsolve)):
                ch1[m, n] = freqtemp
                las.change_freq(3, freqtemp)
        elif mold == m:
            """Afterwards it proceeds with a solve using the faster method. Unfortunately, this
            method may skip to the other solution if conditions are (un)favorable.   (Un)favorable conditions
            include heavy oscillations and bad initial guess for the isclose value."""
            freqsolve, amt = pc.phasematch.solve_frequency(
                samp1, las, 2, 3, isclose=True, amt=amount
            )
            if np.any(list(freqsolve)):
                ch1[m, n] = list(freqsolve)[0]
                las.change_freq(3, list(freqsolve)[0])
            else:
                ch1[m, n] = float("nan")
        else:
            """This final step is testing whether it is better to use the original solve upon a
            new scanline or to stick with the recent solve.  Currently in place is to roll back to
            the original solve.  It then updates angletemp for the next scanline."""
            las.change_freq(3, freqtemp)
            freqsolve, amt = pc.phasematch.solve_frequency(
                samp1, las, 2, 3, isclose=True, amt=amount
            )
            mold = m
            if np.any(list(freqsolve)):
                ch1[m, n] = list(freqsolve)[0]
                freqtemp = list(freqsolve)[0]
                las.change_freq(3, list(freqsolve)[0])
            else:
                ch1[m, n] = float("nan")

for m in range(len(var1a)):
    for n in range(len(var2a)):
        las.change_freq(2, var1a[m])
        las.change_freq(1, var2a[n])
        las.change_freq(3, ch1[m, n])
        Mlist, Mphase, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
        test1[m, n] = -np.log10(Mlist[1])
        # test1[m, n] = (Mlist[1])


data = wt.Data(name="freq3 solves")

data.create_variable(name="w1", units="wn", values=var1)
data.create_variable(name="w2", units="wn", values=var2)
data.create_channel(name="w3freq", values=ch1)
data.create_channel(name="test", values=test1)
data.transform("w2", "w1")

data.channels[0].null = np.min(ch1)

wt.artists.quick2D(data, channel=0)
plt.show()

wt.artists.quick2D(data, channel=1)
plt.show()

pass
