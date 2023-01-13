import numpy as np

# import _layer as Layer
import json
from types import SimpleNamespace


class IsoSample:
    def __init__(self, description=None):
        if description is None:
            self.description = ""
        else:
            self.description = description
        self.layers = list()

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key):
        return setattr(self, key)

    def create_layer(self, csvfilelist, molfraclist, wspacing=1.0, thickness=0.01, label=""):
        """create layer from a list of csvfiles, each file must contain three columns (L-R):  freq (cm-1), absorption coeff (cm-1), n
        Freqs should be in increasing order.

        Parameters:
        -----------
        csvfilelist:  list (str)
            list of string paths to tab-delimited spreadsheet files
        molfraclist:  list (float)
            mole fractions of each component, which should sum to 1.000 but is not checked as such
        thickness: float
           thickness of layer in cm
        label: str
           description of layer

        """
        wbeg = list()
        wend = list()
        wdatalist = list()
        adatalist = list()
        ndatalist = list()

        for i in range(len(csvfilelist)):
            csvfile = csvfilelist[i]
            data = np.loadtxt(csvfile)
            w_i = np.asarray(data[:, 0], dtype=float)
            a_i = np.asarray(data[:, 1], dtype=float)
            n_i = np.asarray(data[:, 2], dtype=float)

            if w_i[0] > w_i[1]:
                return IndexError(f"freqs for list object {i} must be increasing order")

            wbeg.append(w_i[0])
            wend.append(w_i[len(w_i) - 1])
            wdatalist.append(w_i)
            adatalist.append(a_i)
            ndatalist.append(n_i)

        wbegt = np.max(wbeg)
        wendt = np.min(wend)
        number = int((wendt - wbegt) / wspacing) + 1
        wvec = np.linspace(wbegt, wendt, number)

        a2array = np.zeros([len(wvec), i + 1], dtype=float)
        n2array = np.zeros([len(wvec), i + 1], dtype=float)

        for n in range(len(wvec)):
            w = wvec[n]

            for i in range(len(csvfilelist)):
                wi = wdatalist[i]
                ai = adatalist[i]
                ni = ndatalist[i]

                acalc1 = np.interp(w, wi, ai)
                ncalc1 = np.interp(w, wi, ni)

                a2array[n, i] = acalc1 * molfraclist[i]
                n2array[n, i] = ncalc1 * molfraclist[i]

        aout = np.sum(a2array, 1)
        nout = np.sum(n2array, 1)

        layer = Layer()
        layer["label"] = label
        layer["thickness"] = thickness
        layer["w_points"] = wvec
        layer["a_points"] = aout
        layer["n_points"] = nout
        self.layers.append(layer)
        return 0

    def load_layer(self, csvfile, thickness, label=""):
        """load a layer from a tab-delimited spreadsheet file
        File must contain three columns (L-R):  freq (cm-1), absorption coeff (cm-1), n
        Freqs should be in increasing order.

        Parameters:
        -----------
        csvfile: path
           path to tab-delimited spreadsheet file
        thickness: float
           thickness of layer in cm
        label: str
           description of layer (str)
        """
        data = np.loadtxt(csvfile)
        wp = np.asarray(data[:, 0], dtype=float)

        if wp[0] > wp[1]:
            return IndexError("freqs must be increasing order")
        layer = Layer()
        layer["label"] = label
        layer["thickness"] = thickness
        layer["w_points"] = wp
        layer["a_points"] = np.asarray(data[:, 1], dtype=float)
        layer["n_points"] = np.asarray(data[:, 2], dtype=float)
        self.layers.append(layer)
        return 0

    def change_layer(self, layernum, csvfile=None, thickness=None, label=None):
        """Replace a layer with the given number as per the csvfile, thickness, and label.

        Parameters:
        -----------
        layernum : int
           layer number to change
        csvfile: path (optional)
           path to tab-delimited spreadsheet file
        thickness: float (optional)
           thickness of layer in cm
        label: str (optional)
           description of layer (str)
        """

        if csvfile is not None:
            data = np.loadtxt(csvfile)
            wp = np.asarray(data[:, 0], dtype=float)

            if wp[0] > wp[1]:
                return IndexError("freqs must be increasing order")

            self.layers[layernum - 1]["w_points"] = wp
            self.layers[layernum - 1]["a_points"] = np.asarray(data[:, 1], dtype=float)
            self.layers[layernum - 1]["n_points"] = np.asarray(data[:, 2], dtype=float)

        if thickness is not None:
            self.layers[layernum - 1]["thickness"] = thickness

        if label is not None:
            self.layers[layernum - 1]["label"] = label

        return 0

    def as_dict(self):
        """dict representation for this IsoSample class."""
        out = dict()
        out["description"] = self.description
        out["layers"] = list()
        for k in self.layers:
            out["layers"].append(k.as_vars())
        return out

    def save(self, file):
        """Save the JSON representation into an open file."""

        class DictEncoder(json.JSONEncoder):
            def default(self, obj):
                if hasattr(obj, "tolist"):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        f = open(file, "w")
        json.dump(self.as_dict(), f, cls=DictEncoder)
        f.close()
        return 0

    def load(self, file):
        """
        Load a file as a string, convert to IsoSample object
        """
        f = open(file, "r")
        str1 = f.read()
        f.close()
        obj = SimpleNamespace(**eval(str1))
        return self.to_IsoSample(obj)

    def to_IsoSample(self, nsp):
        self.__init__()
        self.description = nsp.description
        for k in nsp.layers:
            m = Layer()
            m["label"] = str(k["label"])
            m["thickness"] = float(k["thickness"])
            m["w_points"] = k["w_points"]
            m["a_points"] = k["a_points"]
            m["n_points"] = k["n_points"]
            self["layers"].append(m)
        return 0


class Layer:
    def __init__(self):
        self = dict()
        self["label"] = ""
        self["thickness"] = 0.00
        self["w_points"] = list()
        self["n_points"] = list()
        self["a_points"] = list()

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, val):
        return setattr(self, key, val)

    def as_vars(self):
        class DictEncoder(json.JSONEncoder):
            def default(self, obj):
                if hasattr(obj, "tolist"):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)

        return vars(self)

    def estimate(self, freq):
        """using the points in the layer, interpolate to estimate the refractive index and
        # absorption coefficient for the specified frequency (cm-1)

        Parameters:
        -----------
        freq: float
           frequency (cm-1)

        Returns:
        ---------
            tuple:  (freq, absorp , ncalc )
                freq : float
                    frequency (cm-1) looped back
                absorp : float
                    absorption coefficient (cm-1)
                ncalc : float
                    real refractive index
        """
        freq1 = float(freq)
        ncalc = np.interp(freq1, self.w_points, self.n_points)
        absorp = np.interp(freq1, self.w_points, self.a_points)

        return freq1, absorp, ncalc

    def suppress_absorbances(self):
        """zero out all absorbance data in the layer."""

        abszero = np.zeros(len(self.a_points))
        self.a_points = abszero

        return 0
