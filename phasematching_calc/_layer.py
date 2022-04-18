import numpy as np
import json

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
        """ using the points in the layer, interpolate to estimate the refractive index and
        # absorption coefficient for the specified frequency (cm-1)

        # Parameters:
        # -----------
        # freq: float
        #    frequency (cm-1) to estimate for
        #
        # Returns:
        # tuple: (freq, absorp, n)
        #       freq : float
        #           loopedback frequency(cm-1)
        #       absorp : float
        #           absorption coefficient (cm-1),
        #       n : float
        #           real refractive index
        """
        freq1=float(freq)
        ncalc=np.interp(freq1,self.w_points,self.n_points)
        absorp=np.interp(freq1,self.w_points,self.a_points)

        return freq1, absorp, ncalc


    def suppress_absorbances(self):
        """ zero out all absorbance data in the layer."""

        abszero=np.zeros(len(self.a_points))
        self.a_points=abszero

        return 0

