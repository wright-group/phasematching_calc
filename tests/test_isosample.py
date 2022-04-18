import phasematching_calc as pc
import os
from config.definitions import ROOT_DIR

"""Test Requires R/W permissions in tests folder."""


filepath=os.path.join(ROOT_DIR, 'tests')

jsonfile=os.path.join(filepath, 'samp1.json')
lay1file=os.path.join(filepath, 'sapphire1.txt')
lay2file=os.path.join(filepath, 'H2O_1.txt')

samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.load_layer(lay1file, 0.02, label="sapphirefw")
samp1.load_layer(lay2file, 0.01, label="water")
samp1.load_layer(lay1file, 0.02, label="sapphirebw")

samp1.save(jsonfile)

savefile=os.path.join(filepath, 'testsum.txt')

samp2=pc.IsoSample.IsoSample()
samp2.load(jsonfile)
assert samp2['description']==desc

