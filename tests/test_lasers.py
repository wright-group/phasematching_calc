import phasematching_calc as pc
import os


"""Test Requires R/W permissions in tests folder."""


def test_laser_save_load():

    filepath = os.path.join(os.getcwd(), "tests")
    jsonfile = os.path.join(filepath, "las.json")

    las = pc.Lasers()
    arr1 = [1500.0, 3000.0, 15000.0]
    las.add_frequencies(arr1)
    arr2 = [0.0, -15.0, 15.0]
    las.add_angles(arr2)
    arr3 = [1, -1, 1]
    las.add_k_coeffs(arr3)
    arr4 = [0, 1, 0]
    las.add_pols(arr4)
    las.change_geometry()
    las.save(jsonfile)

    las2 = pc.Lasers()
    las2.load(jsonfile)

    assert isinstance(las2, pc.Lasers)


if __name__ == "__main__":
    test_laser_save_load()
