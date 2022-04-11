import phasematching_calc as pc
import os
from config.definitions import ROOT_DIR

"""Test Requires R/W permissions in tests folder."""

filepath=os.path.join(ROOT_DIR, 'tests')

jsonfile=os.path.join(filepath, 'las.json')

las=pc.Lasers.Lasers()
arr1=[1500.0,3000.0,15000.0]
las.add_frequencies(arr1)
arr2=[0.0,-15.0, 15.0]
las.add_angles(arr2)
arr3=[1,-1,1]
las.add_k_coeffs(arr3)
arr4=[0,1,0]
las.add_pols(arr4)
las.change_geometry()
las.save(jsonfile)

las2=pc.Lasers.Lasers()
las2.load(jsonfile)

assert (isinstance(las2, pc.Lasers.Lasers))