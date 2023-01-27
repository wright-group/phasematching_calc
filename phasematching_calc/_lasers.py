import numpy as np
import json

__all__ = ["Lasers"]


class Lasers:
    def __init__(self, description=None):
        if description is None:
            self.description = ""
        else:
            self.description = description
        self.frequencies = list()
        self.anglesairdeg = list()
        self.k_coeffs = list()
        self.polarizations = list()
        self.geometry = "boxcars"
        self.supportedgeometrylist = {"boxcars", "planar"}
        self.xmask = [1, 0, 0, 1]
        self.ymask = [0, 1, 1, 0]
        return

    def as_dict(self):
        """Dictionary representation for this Lasers class."""
        out = {}
        out["description"] = self.description
        out["frequencies"] = self.frequencies
        out["anglesairdeg"] = self.anglesairdeg
        out["kcoeffs"] = self.k_coeffs
        out["polarizations"] = self.polarizations
        out["geometry"] = self.geometry
        return out

    def add_frequencies(self, ar1):
        """Add an array of frequencies(cm-1) associated with the laser frequency.
        Value:  float.  Must be of length equal to input frequencies. Current geometries allow only three distinct inputs
        indicative that the array must be of length 3.

        Parameter
        --------
        ar1 : list(float)
            array of frequencies currently only length 3
        """
        ar2 = np.asarray(ar1, dtype=float)
        self.frequencies = ar2
        if len(ar2) != 3:
            return ValueError(
                "geometries require three lasers only(set k coeff to zero if not used)"
            )
        if len(self.anglesairdeg) == 3:
            return self.calculate_cartesian_angles()
        else:
            return 0

    def add_angles(self, ar1):
        """Add an array of angles(degrees) associated with the laser frequency and its geometry.
        These angles correspond to the angles in air
        That the inputs make being focused into the sample.

        Parameter
        --------
        ar1 : list(float)
            arrar of angles (degrees) in air currently only length 3
        """

        ar2 = np.asarray(ar1, dtype=float)
        self.anglesairdeg = ar2
        if len(ar2) != 3:
            return ValueError(
                "geometries require three lasers only(set k coeff to zero if not used)"
            )
        if len(self.frequencies) == 3:
            return self.calculate_cartesian_angles()
        else:
            return 0

    def add_k_coeffs(self, ar1):
        """Add an array of coefficients equivalent to the number and sign of wavevector associated with that input laser frequency.
        Value:  must be int, and can be zero.  Must be of length equal to input frequencies.

        Parameter
        --------
        ar1 : list(int)
            coefficient for the wavemixing for the inputs, currently only length 3
        """
        if type(ar1[0]) != int:
            return ValueError("List elements are not integers")
        ar2 = np.asarray(ar1, dtype=int)
        if len(ar2) != 3:
            return ValueError(
                "geometries require three lasers only(set k coeff to zero if not used)"
            )
        self.k_coeffs = ar2
        return 0

    def add_pols(self, ar1):
        """Add an array of polarizations with length equal to the input frequencies.
        Value:  1 == Vertical,   !1 == Horizontal.   No other polarizations supported, and must be int.

        Parameter
        --------
        ar1 : list(int)
            coefficient for polarizations, currently only length 3
        """
        if type(ar1[0]) != int:
            return ValueError("Polarizations not supported beyond int")
        ar2 = np.asarray(ar1, dtype=int)
        if len(ar2) != 3:
            return ValueError(
                "geometries require three lasers only(set k coeff to zero if not used)"
            )
        self.polarizations = ar2
        return 0

    def change_freq(self, pos, newval):
        """changefreq(self,pos,newval)
        Change the laser frequency of one of the inputs to a new value (cm-1)

        Parameters
        --------
        pos : int
            input number
        newval : float
            frequency to change to (cm-1)
        """
        val = float(newval)
        if val < 0.00:
            return ValueError("must be real greater than zero")
        lenfreq = len(self.frequencies)
        if pos > lenfreq:
            return IndexError("position greater than list length")
        self.frequencies[pos - 1] = val
        return 0

    def change_angle(self, pos, newval):
        """changeangle(self,pos,newval)
        Change the angle of one of the inputs to a new value (degrees)

        Parameters
        --------
        pos : int
            input number
        newval : float
            angle in air (degrees) to change to
        """

        lenfreq = len(self.frequencies)
        if pos > lenfreq:
            return IndexError("position greater than list length")
        self.anglesairdeg[pos - 1] = float(newval)
        return self.calculate_cartesian_angles()

    def change_geometry(self, newval="boxcars"):
        """Change Lasers geometry to one supported by the list.  See alo Lasers.supportedgeometrylist.
        "boxcars":                           "planar":

                  O 2
                  |                                                                              |  y
         4        |         1                                                                    |
         O--------+---------O                ---------+----O-----O--O-------O                    +------ x
                  |                                        2     3  4       1
                  |
                  O 3

        where "+" is the center of focus of the beams, and 4 is the output location.

        Parameters
        --------
        newval : str
            geometry to change to
        """

        if newval in self.supportedgeometrylist:
            self.geometry = newval
            self.calculate_cartesian_angles()
            return self.calculate_masks()
        else:
            return ValueError(
                "unsupported geometry (see supportdgeometrylist for full list of supported geometries)"
            )

    def calculate_cartesian_angles(self):
        """Based on the geometry, convert the angles into x and y angles, in radians."""

        num = len(self.frequencies)
        self.anglesxrad = np.zeros(num)
        self.anglesyrad = np.zeros(num)
        angleairs = (
            np.array(self.anglesairdeg, dtype=float) / 180.000 * np.pi
        )  # converts to radians
        if num == 3:
            if self.geometry == "boxcars":
                self.anglesxrad[0] = angleairs[0]
                self.anglesyrad[1] = angleairs[1]
                self.anglesyrad[2] = angleairs[2]  # NOTE
            elif self.geometry == "planar":
                self.anglesxrad[0] = angleairs[0]
                self.anglesxrad[1] = angleairs[1]
                self.anglesxrad[2] = angleairs[2]
            else:
                return ValueError(
                    "unsupported geometry (see supportdgeometrylist for full list of supported geometries)"
                )
        else:
            return ValueError(
                "geometries require three lasers only(set k coeff to zero if not used"
            )
        return 0

    def calculate_masks(self):
        """Based on the supported geometry list, generate a 0,1 binary type array of x and ys. If 0, that element is
        not to be used in related phasematching calculations.  Binary values are based on whether the inputs
        need that coordinate specified.  For example, on boxcars, input 2 would have a binary 0 for x and binary 1 for y."""
        geom = self.geometry
        if geom == "boxcars":
            self.xmask = [1, 0, 0, 1]
            self.ymask = [0, 1, 1, 0]
        elif geom == "planar":
            self.xmask = [1, 1, 1, 1]
            self.ymask = [0, 0, 0, 0]
        return 0

    def save(self, file):
        """Save the JSON representation into an open file."""
        f = open(file, "w")

        class NdarrayEncoder(json.JSONEncoder):
            def default(self, obj):
                if hasattr(obj, "tolist"):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        json.dump(self.as_dict(), f, cls=NdarrayEncoder)
        f.close()
        return 0

    def load(self, file):
        """Load the JSON representation into a Lasers object."""

        class NdarrayDecoder(json.JSONDecoder):
            def default(self, obj):
                if hasattr(obj, "tolist"):
                    return obj.tolist()
                return json.JSONDecoder.default(self, obj)

        f = open(file, "r")
        ar1 = json.load(f, cls=NdarrayDecoder)
        f.close()
        self.frequencies = ar1["frequencies"]
        self.anglesairdeg = ar1["anglesairdeg"]
        self.k_coeffs = ar1["kcoeffs"]
        self.polarizations = ar1["polarizations"]
        self.geometry = ar1["geometry"]
        self.calculate_cartesian_angles()
        return self.calculate_masks()
