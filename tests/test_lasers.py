import phasematching_calc as pc

import os
from config.definitions import ROOT_DIR
filepath=os.path.join(ROOT_DIR, 'tests')

jsonfile=os.path.join(filepath, 'las.json')

las=pc._lasers.Lasers()
arr1=[1500.0,3000.0,15000.0]
las.addfrequencies(arr1)
arr2=[0.0,-15.0, 15.0]
las.addangles(arr2)
arr3=[1,-1,1]
las.addkcoeffs(arr3)
arr4=[0,1,0]
las.addpolarizations(arr4)
las.changegeometry()
las.save(jsonfile)

las2=pc._lasers.Lasers()
las2.load(jsonfile)

assert (isinstance(las2, pc._lasers.Lasers))