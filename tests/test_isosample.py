import phasematching_calc as pc

import os
from config.definitions import ROOT_DIR
filepath=os.path.join(ROOT_DIR, 'tests')

jsonfile=os.path.join(filepath, 'samp1.json')
lay1file=os.path.join(filepath, 'sapphire1.txt')
lay2file=os.path.join(filepath, 'H2O_1.txt')

samp1=pc._isosample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.loadlayer(lay1file, 0.02, label="sapphirefw")
samp1.loadlayer(lay2file, 0.01, label="water")
samp1.loadlayer(lay1file, 0.02, label="sapphirebw")

samp1.save(jsonfile)

samp2=pc._isosample.IsoSample()
samp2.load(jsonfile)
assert samp2['description']==desc

