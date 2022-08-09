import phasematching_calc as pc
import WrightTools as wt
import matplotlib.pyplot as plt
import numpy as np
import os
from config.definitions import ROOT_DIR
from sympy import S, FiniteSet, Interval, oo


def test_pm_mcalcs():
    # filepath = os.path.join(ROOT_DIR, "tests")
    filepath = os.path.join(os.getcwd(), "tests")
    lay1file = os.path.join(filepath, "CaF2_Malitson.txt")
    lay2file = os.path.join(filepath, "H2O_1.txt")

    tkcaf2 = 0.02  # cm
    tkwater = 0.01  # cm

    # generation of a IsoSample
    samp1 = pc.IsoSample.IsoSample()
    desc = "FWM cell"
    samp1.description = desc
    samp1.load_layer(lay1file, tkcaf2, label="caf2fw")
    samp1.load_layer(lay2file, tkwater, label="water")
    samp1.load_layer(lay1file, tkcaf2, label="caf2bw")

    # generation of a Lasers object.
    las = pc.Lasers.Lasers()
    arr1 = [2700.0, 1800.0, 1800.0]
    las.add_frequencies(arr1)
    arr2 = [10.0, 10.0, 10.0]
    las.add_angles(arr2)
    arr3 = [1, -1, 1]
    las.add_k_coeffs(arr3)
    arr4 = [1, 1, 1]
    las.add_pols(arr4)
    las.change_geometry("boxcars")

    # Other test scripts show how the above objects can be saved and loaded.
    # single point Mcalc check
    Mlist, Mdeltalist, tklist, Tdict = pc.phasematch.m_calc(samp1, las)
    testoutput = np.abs(Mlist[1])
    assert np.isclose(testoutput, 0.029799686)

    # angle estimation for laser 1 in layer 2 with frequency 2600 cm-1
    # using boxcars..test for all angles
    angle, amt = pc.phasematch.solve_angle(samp1, las, 2, 1, 2600)
    assert np.isclose(
        float(angle.end), float(77.9637423)
    )  # angle is an Interval, Interval.end is a Sympy function

    # frequency estimate for a laser 1 in layer 2 to allow for phasematching
    # with the angle in air shown in the Las object above...test for all reals >0
    freq = pc.phasematch.solve_frequency(samp1, las, 2, 1, 10)
    out = freq
    assert out == Interval(
        0, oo
    )  # this may be updated if critical angles will be calculated at high freq in future

    ####### BEGIN SECTION that can be improved next round ##########
    Alist_in, Alist_out = pc.phasematch.calculate_absorbances(samp1, las)
    assert np.isclose(Alist_in[1][0], 1.180156)
    assert np.isclose(Alist_in[1][1], 2.973818)

    tlist_in, tlist_out = pc.phasematch.calculate_ts(samp1, las)
    assert np.isclose(tlist_in[1][0], 1413.241232)
    assert np.isclose(tlist_out[1], 1406.194641)

    Mlistout = pc.phasematch.apply_absorbances(Mlist, Alist_in, Alist_out)
    assert np.isclose(Mlistout[2], 5.5515857e-15)

    Mlistout2 = pc.phasematch.apply_trans(Mlist, Tdict)
    assert np.isclose(Mlistout2[2], 0.97078472)
    ####### END SECTION that can be improved next round ##########


def test_emptyset():
    # new lasers object:  unphasematchable geometry TSF
    # filepath = os.path.join(ROOT_DIR, "tests")
    filepath = os.path.join(os.getcwd(), "tests")

    lay1file = os.path.join(filepath, "CaF2_Malitson.txt")
    lay2file = os.path.join(filepath, "H2O_1.txt")

    tkcaf2 = 0.02  # cm
    tkwater = 0.01  # cm

    # generation of a IsoSample
    samp1 = pc.IsoSample.IsoSample()
    desc = "FWM cell"
    samp1.description = desc
    samp1.load_layer(lay1file, tkcaf2, label="caf2fw")
    samp1.load_layer(lay2file, tkwater, label="water")
    samp1.load_layer(lay1file, tkcaf2, label="caf2bw")

    las2 = pc.Lasers.Lasers()
    arr1 = [6200.0, 3000.0, 25000.0]
    las2.add_frequencies(arr1)
    arr2 = [15.0, 5.0, 10.0]
    las2.add_angles(arr2)
    arr3 = [1, 1, 1]
    las2.add_k_coeffs(arr3)
    arr4 = [1, 1, 1]
    las2.add_pols(arr4)
    las2.change_geometry("boxcars")

    angle, amt = pc.phasematch.solve_angle(samp1, las2, 2, 1, 6000)
    out = angle
    assert out == S.EmptySet


if __name__ == "__main__":
    test_pm_mcalcs()
    test_emptyset()
