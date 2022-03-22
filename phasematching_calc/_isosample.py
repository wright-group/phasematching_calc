import numpy as np
#import _layer as Layer
import json
from types import SimpleNamespace

class IsoSample():
    def __init__(self):
        self.description=""
        self.layers=list()

    def __getitem__(self,key):
        return getattr(self,key)

    def __setitem__(self,key):
        return setattr(self,key)

    def loadlayer(self, csvfile, thickness, label=""):
        # load a layer from a tab-delimited spreadsheet file
        # File must contain three columns (L-R):  freq (cm-1), absorption coeff (cm-1), n
        # Freqs should be in increasing order.

        # Parameters:
        # -----------
        # csvfile: spreadsheet file
        # thickness:  thickness of layer in cm
        # label:  description of layer (str)
        data=np.loadtxt(csvfile)
        wp=np.asarray(data[:,0], dtype=float)

        if (wp[0] > wp[1]):
            return IndexError('freqs must be increasing order')
        layer=Layer()
        layer['label']=label
        layer['thickness']=thickness
        layer['w_points']=wp
        layer['a_points']=np.asarray(data[:,1], dtype=float)
        layer['n_points']=np.asarray(data[:,2], dtype=float)
        self.layers.append(layer)
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

        f=open(file, 'w')
        json.dump(self.as_dict(), f, cls=DictEncoder)
        f.close()
        return 0


    def load(self, file):
        '''
        Load a file as a string, convert to IsoSample object
        '''
        f=open(file,'r')
        str1=f.read()
        f.close()
        obj=SimpleNamespace(**eval(str1))
        return self.to_IsoSample(obj)

        

    def to_IsoSample(self, nsp):
        self.__init__()
        self.description=nsp.description
        for k in nsp.layers:
            m=Layer()
            m['label']=str(k['label'])
            m['thickness']=float(k['thickness'])
            m['w_points']=k['w_points']
            m['a_points']=k['a_points']
            m['n_points']=k['n_points']    
            self["layers"].append(m)
        return 0

class Layer:
    def __init__(self):
        self=dict()
        self['label']=""
        self['thickness']=0.00
        self['w_points']=list()
        self['n_points']=list()
        self['a_points']=list()


    def __getitem__(self,key):
        return getattr(self,key)


    def __setitem__(self,key,val):
        return setattr(self,key,val)


    def as_vars(self):
        class DictEncoder(json.JSONEncoder):
            def default(self, obj):
                if hasattr(obj, "tolist"):
                    return obj.tolist()
                return json.JSONEncoder.default(self, obj)
        return vars(self)


    def estimate(self, freq):
        # using the points in the layer, interpolate to estimate the refractive index and 
        # absorption coefficient for the specified frequency (cm-1)
        
        # Parameters:
        # -----------
        # freq: frequency (cm-1)
        #
        # Returns:
        # tuple: {refractive index, absorption coefficient (cm-1)}
        freq1=float(freq)
        ncalc=np.interp(freq1,self.w_points,self.n_points)
        absorp=np.interp(freq1,self.w_points,self.a_points)

        return freq1, absorp, ncalc