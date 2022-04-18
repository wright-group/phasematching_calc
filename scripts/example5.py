import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import *


filepath = os.path.join(ROOT_DIR, "tests")

lay3file = os.path.join(filepath, "CaF2_Malitson.txt")
lay4file = os.path.join(filepath, "CH3CN_paste_1.txt")
lay5file = os.path.join(filepath, "CaF2_Malitson.txt")

tkcaf2 = 0.02  # cm
tkacn = 0.1  # cm


# generation of a IsoSample
samp1 = pc.IsoSample.IsoSample()
desc = "FWM cell"
samp1.description = desc
samp1.load_layer(lay5file, tkcaf2, label="caf2fw")
samp1.load_layer(lay4file, tkacn, label="ACN")
samp1.load_layer(lay3file, tkcaf2, label="caf2bw")

# new Lasers object
las4 = pc.Lasers.Lasers()
arr1 = [3150.0, 2200.0, 20000.0]
las4.add_frequencies(arr1)
arr2 = [5.0, 10.0, 0.0]
las4.add_angles(arr2)
arr3 = [1, -1, 1]
las4.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las4.add_pols(arr4)
las4.change_geometry("planar")


tin, tout = pc.phasematch.calculate_ts(samp1, las4)
print(tin, tout)

for m in range(len(tin)):
    if m == 0:
        pass
    else:
        for i in range(len(tin[m])):
            tin[m][i] = tin[m][i] - tin[m - 1][i]

for i in range(len(tout)):
    if i == 0:
        pass
    else:
        tout[i] = tout[i] - tout[i - 1]

print(tin, tout)
tlist = list()
x1 = list()
x2 = list()
x3 = list()
x4 = list()
y1 = list()
y2 = list()
y3 = list()
y4 = list()

for m in range(len(tin)):
    tinvec = list(tin[m])
    tinvec.append(tout[m])
    avg = np.mean(tinvec)
    tinvec = np.asarray(tinvec - avg)
    for i in range(len(tinvec)):
        if i == 0:
            x1.append(m + 1)
            y1.append(tinvec[i])
        elif i == 1:
            x2.append(m + 1)
            y2.append(tinvec[i])
        elif i == 2:
            x3.append(m + 1)
            y3.append(tinvec[i])
        elif i == 3:
            x4.append(m + 1)
            y4.append(tinvec[i])
        else:
            pass

plt.rcParams["figure.autolayout"] = True
plt.xlim(0, 5)
plt.ylim(-60.0, 30.0)
plt.grid()
# plt.scatter(x1,y1,marker="o",markersize=10, markeredgecolor="red", markerfacecolor="red")
# plt.scatter(x2,y2,marker="o",markersize=10, markeredgecolor="green", markerfacecolor="green")
# plt.scatter(x3,y4,marker="o",markersize=10, markeredgecolor="blue", markerfacecolor="blue")
# plt.scatter(x4,y4,marker="o",markersize=10, markeredgecolor="black", markerfacecolor="black")

print(x1)
print(y1)

print(x2)
print(y2)

xn1 = x1
yn1 = y1
plt.scatter(xn1, yn1, c="red")

xn1 = x2
yn1 = y2
plt.scatter(xn1, yn1, c="green")

xn1 = x3
yn1 = y3
plt.scatter(xn1, yn1, c="blue")

xn1 = x4
yn1 = y4
plt.scatter(xn1, yn1, c="black")


plt.show()

pass
