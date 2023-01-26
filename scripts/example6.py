import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from sympy import *


filepath = os.path.join(os.getcwd(), "tests")
lay3file = os.path.join(filepath, "CaF2_Malitson.txt")
lay4file = os.path.join(filepath, "CH3CN_paste_1.txt")

tkcaf2 = 0.02  # cm
tkacn = 0.01  # cm

# samp1 = pc.IsoSample.IsoSample()
samp1 = pc.IsoSample()
desc = "FWM cell with fw caf2, sample acn, and bw caf2"
samp1.description = desc
samp1.load_layer(lay3file, tkcaf2, label="caf2fw")
samp1.load_layer(lay4file, tkacn, label="acn")
samp1.load_layer(lay3file, tkcaf2, label="caf2bw")

# las4 = pc.Lasers.Lasers()
las4 = pc.Lasers()
arr1 = [3150.0, 2200.0, 17200.0]
las4.add_frequencies(arr1)
arr2 = [6.0, -13.20, 0.0]
las4.add_angles(arr2)
arr3 = [1, -1, 1]
las4.add_k_coeffs(arr3)
arr4 = [1, 1, 1]
las4.add_pols(arr4)
las4.change_geometry("planar")

angl1, amt = pc.phasematch.solve_angle(samp1, las4, 2, 2)
out = list(angl1)
print(out)

freq, amt = pc.phasematch.solve_frequency(samp1, las4, 2, 3, 10)
out = list(freq)
print(out)

las4.change_freq(3, out[0])

las4.change_freq(2, 2150.0)
angle, amt = pc.phasematch.solve_frequency(samp1, las4, 2, 3, 10)
out2 = list(angle)
print(out2)

las4.change_freq(3, out[0])
angle, amt = pc.phasematch.solve_angle(samp1, las4, 2, 2, isclose=False)
out3 = list(angle)
print(out3)
