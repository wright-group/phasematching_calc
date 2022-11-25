import phasematching_calc as pc
import os
from config.definitions import ROOT_DIR

"""Test Requires R/W permissions in tests folder."""


def test_save_load():
    # filepath = os.path.join(ROOT_DIR, "tests")
    filepath = os.path.join(os.getcwd(), "tests")
    jsonfile = os.path.join(filepath, "samp1.json")
    lay1file = os.path.join(filepath, "sapphire1.txt")
    lay2file = os.path.join(filepath, "H2O_1.txt")

    samp1 = pc.IsoSample.IsoSample()
    desc = "FWM cell"
    samp1.description = desc
    samp1.load_layer(lay1file, 0.02, label="sapphirefw")
    samp1.load_layer(lay2file, 0.01, label="water")
    samp1.load_layer(lay1file, 0.02, label="sapphirebw")

    samp1.save(jsonfile)

    savefile = os.path.join(filepath, "testsum.txt")

    samp2 = pc.IsoSample.IsoSample()
    samp2.load(jsonfile)
    assert samp2["description"] == desc


def test_create_layer():
    filepath = os.path.join(os.getcwd(), "tests")
    lay1file = os.path.join(filepath, "test_file1.txt")
    lay2file = os.path.join(filepath, "test_file2.txt")

    laylist = list()
    laylist.append(lay1file)
    laylist.append(lay2file)
    molfraclist = [0.75, 0.25]

    samp1 = pc.IsoSample.IsoSample()
    desc = "test sample"
    samp1.description = desc
    samp1.create_layer(laylist, molfraclist, wspacing=1.0, thickness=0.01, label="test layer")

    lay = samp1.layers[0]

    layw = lay["w_points"]
    laya = lay["a_points"]
    layn = lay["n_points"]
    thick = lay["thickness"]
    labl = lay["label"]

    assert samp1["description"] == desc

    assert layw[0] == 600.5
    assert len(layw) == int(8961)
    assert layw[8960] == 9560.5

    assert len(layw) == len(laya)
    assert len(layw) == len(layn)
    assert laya[0] == 3.5
    assert layn[0] == 1.75

    assert thick == 0.01
    assert labl == "test layer"


if __name__ == "__main__":
    test_save_load()
    test_create_layer()
