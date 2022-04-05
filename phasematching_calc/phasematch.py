from distutils.log import error
import numpy as np
from ._lasers import Lasers
from ._isosample import IsoSample
from ._isosample import Layer
from sympy import S, FiniteSet, Interval, oo

Iso=IsoSample()
Las=Lasers()

def _guessoutputpol(polsvec):
    """  Guesses output polarization based on isotropic averaging, given the input polarizations and that they are linear and
    aligned with either x or y but not a factor of both."""
    temppolsvec=np.zeros(len(polsvec))
    for i in range(len(polsvec)):
        if (polsvec[i] != 1):
            temppolsvec[i]=0
        else:
            temppolsvec[i]=1
    s1=np.sum(temppolsvec)
    if (s1==3):
        return 1
    elif (s1==2):
        return 0
    elif (s1==1):
        return 1
    elif (s1==0):
        return 0
    else:
        return ValueError("Polarizations list not supported for current estimated output polarization")
    

def _calculatetrans(narr,anglexarr,angleyarr,polsvec,noutvec,angleoutxvec,angleoutyvec,polout):
    '''
    Uses Fresnel's equations to calculate transmission coefficients between layers in a sample,
    as well as after the last layer

    Assumes magnetic susceptibility of materials involved in layers to be zero.
    
    Parameters
    ---------
    narr =     2D array of n where the 1st D is layer and 2nd are the input freqs
    anglexarr = "  of angles in x (radians) where 1st D is first input angles, THEN number of layers 
    angleyarr = "  of angles in y (radians) "
    noutvec  1D array of output n for all layers
    angleoutxvec= " of output angle in x for all layers AND additional output angle (assumed air)
    angleoutyvec= " of output angle in y for all layers "
    polsvec = array of polarizations for input lasers (int)...see list below for polout)
    polout= output polarization (int) (based on tensor averaging, user must assume knowledge of
    isotropic tensors as calculation does not automatically return errors for an output polarization
    not expected to have signal due to isotropic averaging)
        1=vertical
        !1 = horizontal 


    Return
    ------
    tuple:  Tinx, Tiny, Toutx, Touty
    Tin= Transmission coefficient 2D array of 1D=frequencies, 2D=layers frequency in that layer according
    to the geometry and polarizations found in the Lasers object. Polarizations must be linear and array
    of the form 1 = Vertical and != 1 = Horizontal 
    Tinx is x portion and Tiny is y portion.   The normal list of geometries presume that each input
    is either x or y but not both.  

    
    Tout= transmission coefficient through each layer for the output as well as the final exit into the 
    air (m+1 layer array)
    Toutx is the x potion and Touty is y portion.
    '''
    
    Tinx=list()
    Tiny=list()
    Toutx=list()
    Touty=list()


    for m in range(len(noutvec)):
        if (m==0):
            n1vec=np.ones(len(polsvec))
        else:
            n1vec=narr[m-1]
        n2vec=narr[m]
        anglex1vec=anglexarr[m]
        angley1vec=angleyarr[m]
        anglex2vec=anglexarr[m+1]
        angley2vec=angleyarr[m+1]

        Tinxtemp=np.zeros(len(polsvec))
        Tinytemp=np.zeros(len(polsvec))

        for i in range(len(polsvec)):
            
            n2=n2vec[i]
            anglex2=anglex2vec[i]
            angley2=angley2vec[i]
            pol=polsvec[i]
            n1=n1vec[i]
            anglex1=anglex1vec[i]
            angley1=angley1vec[i]

            if pol==1:
                Rx=(n1*np.cos(anglex1)-n2*np.cos(anglex2))**2/(n1*np.cos(anglex1)+n2*np.cos(anglex2))**2
                Ry=(n1*np.cos(angley2)-n2*np.cos(angley1))**2/(n1*np.cos(angley2)+n2*np.cos(angley1))**2
                Tx=1-Rx
                Ty=1-Ry
            else:
                Rx=(n1*np.cos(anglex2)-n2*np.cos(anglex1))**2/(n1*np.cos(anglex2)+n2*np.cos(anglex1))**2
                Ry=(n1*np.cos(angley1)-n2*np.cos(angley2))**2/(n1*np.cos(angley1)+n2*np.cos(angley2))**2
                Tx=1-Rx
                Ty=1-Ry

            Tinxtemp[i]=Tx
            Tinytemp[i]=Ty
        
        Tinx.append(Tinxtemp)
        Tiny.append(Tinytemp)

    
    for m in range(len(noutvec)):
            n1=noutvec[m]
            anglex2=angleoutxvec[m+1]
            angley2=angleoutyvec[m+1]
            anglex1=angleoutxvec[m]
            angley1=angleoutyvec[m]

            if (m==(len(noutvec)-1)):
                n2=1.000
            else:
                n2=noutvec[m+1]
            
            if polout==1:
                Rx=(n1*np.cos(anglex1)-n2*np.cos(anglex2))**2/(n1*np.cos(anglex1)+n2*np.cos(anglex2))**2
                Ry=(n1*np.cos(angley2)-n2*np.cos(angley1))**2/(n1*np.cos(angley2)+n2*np.cos(angley1))**2
                Tx=1-Rx
                Ty=1-Ry
            else:
                Rx=(n1*np.cos(anglex2)-n2*np.cos(anglex1))**2/(n1*np.cos(anglex2)+n2*np.cos(anglex1))**2
                Ry=(n1*np.cos(angley1)-n2*np.cos(angley2))**2/(n1*np.cos(angley1)+n2*np.cos(angley2))**2
                Tx=1-Rx
                Ty=1-Ry

            Toutx.append(Tx)
            Touty.append(Ty)

    return Tinx, Tiny, Toutx, Touty


def _calculatecriticalangle(Iso, Las, layernum, freqnum, frequency=None):
    """  Determines the critical angle for the layernum-1: layernum boundary.  Replaces value at freqnum 
    by frequency if not None.

    Return
    -----
    n,crit float of refractive index, critical angle (radians) for the layernum-1:layernum boundary,
    or numpy.pi/2 if no critical angle found.
    """
 
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")
    if (layernum < 1):
        return ValueError("layernum cannot be less than 1")
    
    if (frequency is not None):
        Las.frequencies[freqnum-1]=frequency

    freq=Las.frequencies[freqnum-1]

    m = layernum-1
    if (layernum==len(Iso['layers'])):
        w,a,nold=Iso["layers"][m].estimate(freq)
        n=float(1.000)
    else:
        w,a,nold=Iso["layers"][m].estimate(freq)
        w,a,n=Iso["layers"][m+1].estimate(freq)
    
    if (n/nold > 1):
        crit=np.pi/2
    else:
        crit=np.arcsin(n/nold)

    return nold, crit


def _calculateallcriticalangles(Iso, Las, freqnum=None, frequency=None):
    """  Calculates all critical angles for all inputs at all layer boundaries.  Or returns
    numpy pi/2 if there is no critical angle at that combination.  Replaces freqnum with 
    frequency if specified.

    Return
    -----
    out= 2D arrray (1st D is layers m, 2nd is laser inputs i with float of the critical angle
    between m and m+1 for each element  or numpy.pi/2 if no critical angle found.
    """
 
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")
    if (frequency < 0.00):
        return ValueError("frequency cannot be less than 0")
    
    if ((frequency is not None) & (freqnum is not None)):
        Las.frequencies[freqnum-1]=frequency

    mt=len(Iso["layers"])
    it=len(Las.frequencies)
    critarr=np.zeros([mt,it])
    narr=critarr

    for m in range(mt):
        for i in range(it):
            n, crit = _calculatecriticalangle(Iso, Las, m+1, i+1)
            critarr[m,i]=crit
            narr[m,i]=n

    return narr, critarr


def _stackcriticalangles(narr, critarr, layernum, freqnum):
    """  Up to a given layernum, finds the smallest of all critical angles of a freqnum and
    returns it.  This is necessary for the SolveAngle result when a large solution of angles are possible.
   
    Return
    -----
    crit = float angle in air (radians) to result in the smallest critangle in critarr up to layernum for freqnum.
    """

    critvec=list()
    nvec=list()

    for m in range(layernum):
        critvec.append(critarr[m][freqnum-1])
        nvec.append(narr[m][freqnum-1])

    crit=np.pi/2

    for m in range(layernum):
        ncurrent=nvec[layernum-m-1]    
        critcurrent=critvec[layernum-m-1]
        crit=np.minimum(crit,critcurrent)

        if (m != layernum):
            nnew=nvec[layernum-m-1]
        else:
            nnew=float(1.000)
        
        critnew=np.arcsin(nnew/ncurrent*np.sin(crit))

        if (np.isnan(critnew)):
            critnew=np.pi/2
        
        crit=critnew

    return crit


def calculateoriginalcritangle(Iso, Las, layernum, freqnum, frequency=None):
    """calculates the maximum original angle in air to result in the critical angle at the
    layernum : layernum+1 interface (or numpy pi/2 if there are no critical angles, or
    the angle in air of the smallest of any critical angles for that freqnum in any previous layers, 
    converted to its original angle in air)

    Return
    ----
    angle:  float (degrees) in air of critical angle
    """

    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")
    if (layernum < 1):
        return ValueError("layernum cannot be less than 1")

    narr,critarr= _calculateallcriticalangles(Iso, Las, freqnum, frequency)
    angle= _stackcriticalangles(narr, critarr, layernum, freqnum)
    angledeg=angle/np.pi*180.00

    return angledeg


def Mcalc(Iso, Las):
    '''
    Calculates the phase mismatching factors for wavemixing for an IsotropicSample, which can then
    be incorporated into plots or simulations (WrightSim)

    Parameters
    ----
    Iso: an IsotropicSample object
    Las: a Lasers object

    Outputs
    ---
    a tuple Mlist, tklist, Tdict consisting of:
    
    Mlist :  a complex array of phasemismatching factors for wavemixing at the output,
    currently only supporting four wave mixing models with supportedgeometries shown in the Laser object.
    
    tklist : the effective thickness of each layer as pertaining to the launched output wave, i.e.
            layer thickness / cosine(angleout)
    Tdict : a dictionary with entries related to transmission coefficients of the lasers and output
    through the sample layers.  NOTE:  This application has not been performed on the entries in Mlist.
    The coefficients rely on linear polarizations defined in the Laser object.  User must perform
    the multiplication.

    Keys in the Tdict are:
    Tdict['Txin']= float(2Darr) Transmission coefficients of the x-coordinated laser inputs (those
    geometries with beams along x=0 would then have 0 values in the entries), each element per layer
    in the second dimension and each laser input in the 1st dimension
    Tdict['Tyin']= " similar for the y-coordinated inputs
    Tdict['Txout']= float (1D) Transmission coefficients of the x-coordinates of the output (those
    geometries with outputs along x=0 would then have 0 values in the entries), each element per layer, with 
    an additional value representing the final transmission coefficient back into air.
    Tdict['Tyout'] = similar for the y-coordinates of the output
    Tdict['launchanglexdeg']= launch angle of output in degrees in air, for xcoordinated output
    Tdict['launchangleydeg']= similar for ycoordinated output
    '''
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")

    numlayers=len(Iso['layers'])
    freqs=Las.frequencies
    kcoeffs=Las.k_coeffs
    pols=Las.polarizations
    anglexrad=Las.anglesxrad
    angleyrad=Las.anglesyrad
    
    numfreqs=len(freqs)
    freqout=float(0.00)
    for m in range(numfreqs):
        freqout=freqout+kcoeffs[m]*freqs[m]
    
    # These are 2D arrays where the 1st D is layer (1st and last are air) and 2nd are the input freqs
    anglex=list()
    angley=list()
    nvec=list()

    # These are 2D arrays where the 1st D is layer (no air layers) and 2nd are the input freqs
    kx=list()
    ky=list()
    kz=list()
    avec=list()
    
    # These are 1D arrays where the D is layer (1st and last are air)
    nout=list()
    angleoutx=list()
    angleouty=list()

    # These are 1D arrays where the D is layer (no air layers)
    tk=list()
    aout=list()
    Mlist=list()
    tklist=list()

    anglex1=anglexrad
    anglex.append(anglex1)

    angley1=angleyrad
    angley.append(angley1)

    for m in range(numlayers):
        anglextemp=np.zeros(numfreqs)
        angleytemp=np.zeros(numfreqs)
        kxtemp=np.zeros(numfreqs)
        kytemp=np.zeros(numfreqs)
        kztemp=np.zeros(numfreqs)
        atemp=np.zeros(numfreqs)
        nvectemp=np.zeros(numfreqs)
        layertemp=Iso['layers'][m]
                    
        wout,aouttemp,nouttemp=layertemp.estimate(freqout)

        angleoutytemp=0.00
        angleoutxtemp=0.00
        koutz=0.00
        kouty=0.00
        koutx=0.00

        for i in range(numfreqs):
            anglex1temp,angley1temp=Angle(Iso,Las,m+1,i+1)
            w,a,n=layertemp.estimate(freqs[i])

            # NOTE: due to specific geometries used so far, it is unnecessary to have
            # r, theta, phi conversions
            #
            if (anglex1temp==0.00):
                anglez=np.pi/2-angley1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                kouty=2*np.pi*n*w*np.sin(angley1temp)*kcoeffs[i]+kouty
                
            else:
                anglez=np.pi/2-anglex1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                koutx=2*np.pi*n*w*np.sin(anglex1temp)*kcoeffs[i]+koutx

            kx1=2*np.pi*n*w*np.sin(anglex1temp)
            ky1=2*np.pi*n*w*np.sin(angley1temp)
            kz1=2*np.pi*n*w*np.sin(anglez)
            anglextemp[i]=anglex1temp
            angleytemp[i]=angley1temp
            kxtemp[i]=kx1
            kytemp[i]=ky1
            kztemp[i]=kz1
            nvectemp[i]=n
            atemp[i]=a

        angleoutxtemp=np.arctan(koutx/koutz)
        angleoutytemp=np.arctan(kouty/koutz)
        anglex.append(anglextemp)
        angley.append(angleytemp)      
        kx.append(kxtemp)
        ky.append(kytemp)
        kz.append(kztemp)
        avec.append(atemp)
        nvec.append(nvectemp)

        tk.append(layertemp['thickness'])
        aout.append(aouttemp)
        angleouty.append(angleoutytemp)
        angleoutx.append(angleoutxtemp)
        nout.append(nouttemp)

    launchanglex=np.arcsin(nouttemp*np.sin(angleoutxtemp))
    launchanglexdeg=launchanglex*180.00/np.pi
    launchangley=np.arcsin(nouttemp*np.sin(angleoutytemp))
    launchangleydeg=launchangley*180.00/np.pi

    angleoutx.append(launchanglex)
    angleouty.append(launchangley)

    polout=_guessoutputpol(pols)
    
    Txin, Tyin, Txout, Tyout = _calculatetrans(nvec,anglex,angley,pols,nout,angleoutx,angleouty,polout)
    Mctemp=1

    for m in range(numlayers):
        avectemp=avec[m]
        kztemp=kz[m]
        kytemp=ky[m]
        kxtemp=kx[m]
        
        tktemp=tk[m]
                
        aouttemp=aout[m]
        angleoutxtemp=angleoutx[m]
        angleoutytemp=angleouty[m]
        nouttemp=nout[m]

        kout=2.00*np.pi*freqout*nouttemp

        if (angleoutxtemp==0):
            tkeff=tktemp/np.cos(angleoutytemp)
            koutx=0
            koutz=kout*np.cos(angleoutytemp)
            koutx=kout*np.sin(angleoutxtemp)
        else:
            tkeff=tktemp/np.cos(angleoutxtemp)
            kouty=0
            koutz=kout*np.cos(angleoutytemp)
            kouty=kout*np.sin(angleoutytemp)

               
        ksumx=ksumy=ksumz=0
        
        for i in range(numfreqs):
            ksumx=kcoeffs[i]*kxtemp[i]+ksumx
            ksumy=kcoeffs[i]*kytemp[i]+ksumy
            ksumz=kcoeffs[i]*kztemp[i]+ksumz
        k4=np.sqrt(ksumx**2+ksumy**2+ksumz**2)
        #dk=np.sqrt((koutx-ksumx)**2+(kouty-ksumy)**2+(koutz-ksumz)**2)
        dk=k4-kout        
        da=0.5*(aouttemp-(np.abs(kcoeffs[0])*avectemp[0]+np.abs(kcoeffs[1])*avectemp[1]+np.abs(kcoeffs[2])*avectemp[2])*tkeff)
        Mc1=np.exp(-0.5*aouttemp*tkeff)
        Mc2=((1-np.exp(da))**2+4*np.exp(da)*(np.sin(dk*tkeff/2))**2)/(da**2+(dk*tkeff)**2)
        
        #Mc2=np.complex(np.cos(dk*tkeff),np.sin(dk*tkeff))*np.exp(-da*tkeff))
        #Mc3=np.complex(-da*tkeff/((dk**2+da**2)*tkeff**2), -tkeff*dk/((dk**2+da**2)*tkeff**2))
        
        #if (i==0):
        #    Mpre=1
        #else:
        #    Mpre=Mc1*Mc2    

        #Mctemp=Mpre*(Mc1)*(Mc2-1)*Mc3
        Mctemp=Mc1*Mc2
        Mlist.append(Mctemp)
        tklist.append(tkeff)

    #angleoutxtemp and angleoutytemp can be converted to launch angles via Snell's Law to calculate
    #launch angle in air afterwards
    launchanglexdeg=np.arcsin(nouttemp*np.sin(angleoutxtemp))*180.00/np.pi
    launchangleydeg=np.arcsin(nouttemp*np.sin(angleoutytemp))*180.00/np.pi

    Tdict=dict()
    Tdict['Txout']=Txout
    Tdict['Tyout']=Tyout
    Tdict['Txin']=Txin
    Tdict['Tyin']=Tyin
    Tdict['launchanglexdeg']=launchanglexdeg
    Tdict['launchangleydeg']=launchangleydeg
    return Mlist, tklist, Tdict


def Angle(Iso,Las,layernum,freqnum, frequency=None):
    '''
    Uses Snell's Law to calculate the angle a laser makes in a specific layer of a sample,
    the laser defined by freqnum and layer by layernum

    Parameters
    ---------
    Iso = IsoSample object
    Las = Lasers object
    layernum = layer number
    freqnum = frequency number (index number)
    (optional) frequency =if not None replaces the frequency at that number with the given value
     prior to calculation

    Output
    ------
    tuple: {anglex,angley} projected angles for the selected frequency in that layer according
    to the geometry found in the Lasers object
    '''
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")
    if (layernum < 1):
        return ValueError("layernum cannot be less than 1")

    lasfreq=Las.frequencies[freqnum-1]
    lasanglex=Las.anglesxrad[freqnum-1]
    lasangley=Las.anglesyrad[freqnum-1]
 
    if (frequency is not None):
        lasfreq=frequency

    for m in range(layernum):
        if (m==0):
            nold=float(1.000)
            aold=float(0.000)
        else:
            nold=n
            aold=a
        w,a,n=Iso["layers"][m].estimate(lasfreq)
        lasanglex=np.arcsin(nold/n*np.sin(lasanglex))
        lasangley=np.arcsin(nold/n*np.sin(lasangley))
    return lasanglex,lasangley


def SolveAngle(Iso,Las,layernum,freqnum, frequency=None):
    '''
    Given an Isotropic Sample, a layer number and frequency number used in a geometry defined in Las object,
    determine the input angle (in air) required for phasematching that frequency. Uses Sympy Set. 
    Returns an empty FiniteSet if a solution cannot be found within an internal convergence iteration series,
    usually meaning there is no solution.
    
    Parameters
    ----
    Iso = The IsotropicSample object
    Las = The Lasers object with only supportedgeometry list capable of solutions
    layernum = layer in which to solve for angle
    freqnum = laser position (defined by geometry)
    frequency = frequency at that laser position

    All other frequencies, k coefficients, refractive indexes are taken from the objects.

    Output
    ----
    Sympy: FiniteSet : theta = (deg) float of original angle in air needed for that PM condition.
           FiniteSet is empty if a solution cannot be found
           Interval(0,90) (deg) if all real angles greater than zero are found, or a more restricted
           interval if restricted by a critical angle before the layernum
    '''
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")

    freqs=Las.frequencies
    kcoeffs=Las.k_coeffs
    
    numfreqs=len(freqs)
    freqout=float(0.00)
    
    flag=int(0)
    
    if (frequency is not None):
        freqs[freqnum-1]=frequency
        if (frequency < 0.00):
            return ValueError("frequency cannot be less than 0")

    for i in range(numfreqs):
        freqout=freqout+kcoeffs[i]*freqs[i]


    for m in range(layernum):
        anglextemp=np.zeros(numfreqs)
        angleytemp=np.zeros(numfreqs)
        kxtemp=np.zeros(numfreqs)
        kytemp=np.zeros(numfreqs)
        kztemp=np.zeros(numfreqs)
        atemp=np.zeros(numfreqs)
        nvectemp=np.zeros(numfreqs)
        layertemp=Iso.layers[m]
        koutz=0.00
        kouty=0.00
        koutx=0.00
        

        for i in range(numfreqs):
            anglex1temp,angley1temp=Angle(Iso,Las,m+1,i+1)
            w,a,n=layertemp.estimate(freqs[i])
            # NOTE: see above
            #
            if (anglex1temp==0.00):
                anglez=np.pi/2.000-angley1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                kouty=2*np.pi*n*w*np.sin(angley1temp)*kcoeffs[i]+kouty
            else:
                anglez=np.pi/2.000-anglex1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                koutx=2*np.pi*n*w*np.sin(anglex1temp)*kcoeffs[i]+koutx
            if (i+1==freqnum):
                factor=0.00
            else:
                factor=1.00
            # this is put in as reminder that we are solving for this variable so the k's for this
            # one have to be set to zero

            kx1=2*np.pi*n*w*np.sin(anglex1temp)*factor
            ky1=2*np.pi*n*w*np.sin(angley1temp)*factor
            kz1=2*np.pi*n*w*np.sin(anglez)*factor
            anglextemp[i]=anglex1temp
            angleytemp[i]=angley1temp
            kxtemp[i]=kx1
            kytemp[i]=ky1
            kztemp[i]=kz1
            nvectemp[i]=n
            atemp[i]=a

        ksumx=ksumy=ksumz=0
        
        for i in range(numfreqs):
            ksumx=kcoeffs[i]*kxtemp[i]+ksumx
            ksumy=kcoeffs[i]*kytemp[i]+ksumy
            ksumz=kcoeffs[i]*kztemp[i]+ksumz
    
        if ((ksumx==ksumz) & (ksumx==0.00)):
            flag=1
            break
        elif ((ksumy==ksumz) & (ksumx==0.00)):
            flag=1
            break


    if (flag==1):
        angledeg=calculateoriginalcritangle(Iso, Las, layernum, freqnum, frequency)
        return Interval(0,angledeg)
    else:
        m = layernum-1
  
        Lastemp=Las
        Isotemp=Iso
        tol=0.0001
        for k in range(layernum):
            Isotemp.layers[k].suppressabs()

        dir=1.00
        amt= 1.00 #deg
        Mtest,tklist,Tdict=Mcalc(Isotemp, Lastemp)

        magMtest=np.abs(Mtest[m]) 
        error1= 1-magMtest
        error2= error1
        iter=70
        angle=Lastemp.anglesairdeg[freqnum-1]
        b=0

        while( error2 > tol):
            b=b+1
            if (error2 > error1):
                dir=(-1.00)*dir
            if (np.abs(error2) < 0.33*np.abs(error1)):
                amt=0.1*amt
            error1=error2
            angle=Lastemp.anglesairdeg[freqnum-1]+amt*dir
            Lastemp.changeangle(freqnum,angle)
            Mtest,tklist,Tdict=Mcalc(Isotemp, Lastemp)
            magMtest=np.abs(Mtest[m]) 
            error2=1-magMtest
            if (b > iter):
                flag=2
                break

        if (flag==2):
            return FiniteSet()
        else:
            return FiniteSet(angle)


def SolveFrequency(Iso, Las, layernum, freqnum, amt=None):
    '''Using the current frequency as first guess, solves for the nearest possible phasematching frequency
    at a fixed angle for that frequencynum in a given layer, using an iterative convergence.  Uses Sympy Set. 
    Returns an empty FiniteSet if a solution cannot be found within an internal convergence iteration series,
    usually meaning there is no solution.
    
    Parameters
    ----
    Iso = The IsotropicSample object
    Las = The Lasers object with only supportedgeometry list capable of solutions
    layernum = layer in which to solve for angle
    freqnum = laser position (defined by geometry)
    amount = amount to change frequency by per convergence step.  Approximates internally if set to None.

    Output
    ----
    Sympy: FiniteSet : frequency = (cm-1) float of frequency needed for that PM condition.
           FiniteSet is empty if a solution cannot be found.
           Interval(0,oo) if all real frequencies greater than zero are found, or a closed
            interval if restricted via some critical angle before the layernum.
    '''
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")
    if (layernum < 1):
        return ValueError("layernum cannot be less than 1")

    freqs=Las.frequencies
    kcoeffs=Las.k_coeffs

    flag=int(0)
    numfreqs=len(freqs)
    freqout=float(0.00)

 
    if amt is None:
        amt = 0.01*freqs[freqnum-1]  # a guess

    for i in range(numfreqs):
        freqout=freqout+kcoeffs[i]*freqs[i]

    for m in range(layernum):
        anglextemp=np.zeros(numfreqs)
        angleytemp=np.zeros(numfreqs)
        kxtemp=np.zeros(numfreqs)
        kytemp=np.zeros(numfreqs)
        kztemp=np.zeros(numfreqs)
        atemp=np.zeros(numfreqs)
        nvectemp=np.zeros(numfreqs)
        layertemp=Iso['layers'][m]
                    
        koutx=0.00
        kouty=0.00
        koutz=0.00

        for i in range(numfreqs):
            anglex1temp,angley1temp=Angle(Iso,Las,m+1,i+1)
            w,a,n=layertemp.estimate(freqs[i])
            # NOTE: see above
            #
            if (i+1 != freqnum):
                if (anglex1temp==0.00):
                    anglez=np.pi/2.000-angley1temp
                    koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                    kouty=2*np.pi*n*w*np.sin(angley1temp)*kcoeffs[i]+kouty
                else:
                    anglez=np.pi/2.000-anglex1temp
                    koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                    koutx=2*np.pi*n*w*np.sin(anglex1temp)*kcoeffs[i]+koutx
            else:
                anglex1temp=0.00
                angley1temp=0.00
                anglez=0.00
            # this is put in as reminder that we are solving for this variable so the k's for this
            # one have to be set to zero

            kx1=2*np.pi*n*w*np.sin(anglex1temp)
            ky1=2*np.pi*n*w*np.sin(angley1temp)
            kz1=2*np.pi*n*w*np.sin(anglez)
            anglextemp[i]=anglex1temp
            angleytemp[i]=angley1temp
            kxtemp[i]=kx1
            kytemp[i]=ky1
            kztemp[i]=kz1
            nvectemp[i]=n
            atemp[i]=a
        
        ksumx=ksumy=ksumz=0

        for i in range(numfreqs):
            ksumx=kcoeffs[i]*kxtemp[i]+ksumx
            ksumy=kcoeffs[i]*kytemp[i]+ksumy
            ksumz=kcoeffs[i]*kztemp[i]+ksumz
    
        if ((ksumx==ksumz) & (ksumx==0.00)):
            flag=1
            break
        elif ((ksumy==ksumz) & (ksumx==0.00)):
            flag=1
            break


    if (flag==1):
        return Interval(0,oo)  #currently does not attempt a walkback of high frequencies
                            # at the given angle to see if it is reflected at a critical angle
    else:
        m = layernum-1
  
        Lastemp=Las
        Isotemp=Iso
        tol=0.0001
        for k in range(layernum):
            Isotemp.layers[k].suppressabs()

        dir=1.00

        Mtest,tklist,Tdict=Mcalc(Isotemp, Lastemp)

        magMtest=np.abs(Mtest[m])
        freq=Lastemp.frequencies[freqnum-1]
        error1= 1-magMtest
        error2= error1
        iter=500
        b=0
        while( error2 > tol):
            b=b+1
            if (error2 > error1):
                dir=(-1.00)*dir
            if (np.abs(error2) < 0.33*np.abs(error1)):
                amt=0.1*amt
            error1=error2
            freq=Lastemp.frequencies[freqnum-1]+amt*dir
            Lastemp.changefreq(freqnum,freq)
            Mtest,tklist,Tdict=Mcalc(Isotemp, Lastemp)
            magMtest=np.abs(Mtest[m])
            error2=1-magMtest
            if (b > iter):
                flag=2
                break

        if (flag==2):
            return FiniteSet()
        else:
            return FiniteSet(freq)


def calculatedeltats(Iso, Las):
    """Calculate the change in delays each pulse makes with respect to the first pulse (first on the Las.frequencies
       list)  as each passes through the IsoSample.

       Uses the refractive indexes, angles of incidence, and thickness per sample, to determine the change in
       timing each pulse makes through a sample relative to the first pulse assuming all are perfectly overlapped
       before entering the IsoSample.

       Is used to confirm that the pulse overlap in thicker samples is satisfactory or not satisfactory enough
       to determine if delay changes are responsible for a limit to the interaction in such a sample.


       Return
       ------
       dt_chart_in:  2D array with frequencies in 1st axis and layernums in second.  Elements are times (fsec)
       in which frequency i  reaches the next layer.     
       dt_chart_out: 1D array with layernums for the output at the given kcoeffs.
       """
    
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")

    numlayers=len(Iso['layers'])
    freqs=Las.frequencies
    kcoeffs=Las.k_coeffs
    pols=Las.polarizations
    anglexrad=Las.anglesxrad
    angleyrad=Las.anglesyrad
    cvac=29979245800  #cm/sec

    numfreqs=len(freqs)
    freqout=float(0.00)
    for m in range(numfreqs):
        freqout=freqout+kcoeffs[m]*freqs[m]
    
    # These are 2D arrays where the 1st D is layer (1st and last are air) and 2nd are the input freqs
    anglex=list()
    angley=list()
    nvec=list()
    
    # These are 1D arrays where the D is layer (1st and last are air)
    nout=list()
    angleoutx=list()
    angleouty=list()

    # These are 1D arrays where the D is layer (no air layers)
    tk=list()

    anglex1=anglexrad
    anglex.append(anglex1)

    angley1=angleyrad
    angley.append(angley1)

    dt_chart_in=np.zeros([numlayers,numfreqs])
    dt_chart_out=np.zeros(numlayers)

    for m in range(numlayers):
        anglextemp=np.zeros(numfreqs)
        angleytemp=np.zeros(numfreqs)
        nvectemp=np.zeros(numfreqs)
        layertemp=Iso['layers'][m]
                    
        wout,aouttemp,nouttemp=layertemp.estimate(freqout)

        angleoutytemp=0.00
        angleoutxtemp=0.00
        koutz=0.00
        kouty=0.00
        koutx=0.00

        for i in range(numfreqs):
            anglex1temp,angley1temp=Angle(Iso,Las,m+1,i+1)
            w,a,n=layertemp.estimate(freqs[i])
            # NOTE: due to specific geometries used so far, it is unnecessary to have
            # r, theta, phi conversions
            #
            if (anglex1temp==0.00):
                anglez=np.pi/2-angley1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                kouty=2*np.pi*n*w*np.sin(angley1temp)*kcoeffs[i]+kouty
              
            else:
                anglez=np.pi/2-anglex1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                koutx=2*np.pi*n*w*np.sin(anglex1temp)*kcoeffs[i]+koutx

            anglextemp[i]=anglex1temp
            angleytemp[i]=angley1temp
            nvectemp[i]=n

        angleoutxtemp=np.arctan(koutx/koutz)
        angleoutytemp=np.arctan(kouty/koutz)
        anglex.append(anglextemp)
        angley.append(angleytemp)      
        nvec.append(nvectemp)
        tk.append(layertemp['thickness'])
        angleouty.append(angleoutytemp)
        angleoutx.append(angleoutxtemp)
        nout.append(nouttemp)

    launchanglex=np.arcsin(nouttemp*np.sin(angleoutxtemp))
    launchangley=np.arcsin(nouttemp*np.sin(angleoutytemp))

    angleoutx.append(launchanglex)
    angleouty.append(launchangley)

    for i in range(numfreqs):
        dttemp=float(0.00)
        for m in range(numlayers):
            anglextemp=anglex[m][i]
            angleytemp=angley[m][i]
            ntemp=nvec[m][i]
            thick=tk[m]
            if (anglextemp==0):
                tkeff=thick/np.cos(angleytemp)
            else:
                tkeff=thick/np.cos(anglextemp)
            dttemp=dttemp+tkeff*ntemp/cvac*1E-15
            dt_chart_in[m][i]=dttemp

    dttemp=float(0.00)
    for m in range(numlayers):
        anglextemp=angleoutx[m][i]
        angleytemp=angleouty[m][i]
        ntemp=nout[m]
        thick=tk[m]
        if (anglextemp==0):
            tkeff=thick/np.cos(angleytemp)
        else:
            tkeff=thick/np.cos(anglextemp)
        dttemp=dttemp+tkeff*ntemp/cvac*1E-15
        dt_chart_out[m]=dttemp

    return dt_chart_in, dt_chart_out


def calculateabsorbances(Iso, Las):
    """Calculate the absorbances each laser (including output) makes per layer in the IsoSample.

       Return (tuple)
       ------
       Alist_in:  2D array with frequencies in 1st axis and layernums in second.  Elements are times (fsec)
       in which frequency i  reaches the next layer.     
       Alist_out: 1D array with layernums for the output at the given kcoeffs.
       """
    
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")

    numlayers=len(Iso['layers'])
    freqs=Las.frequencies
    kcoeffs=Las.k_coeffs
    pols=Las.polarizations
    anglexrad=Las.anglesxrad
    angleyrad=Las.anglesyrad

    numfreqs=len(freqs)
    freqout=float(0.00)
    for m in range(numfreqs):
        freqout=freqout+kcoeffs[m]*freqs[m]
    
    # These are 2D arrays where the 1st D is layer (1st and last are air) and 2nd are the input freqs
    anglex=list()
    angley=list()
    nvec=list()
    avec=list()
    
    # These are 1D arrays where the D is layer (1st and last are air)
    nout=list()
    aout=list()
    angleoutx=list()
    angleouty=list()

    # These are 1D arrays where the D is layer (no air layers)
    tk=list()

    anglex1=anglexrad
    anglex.append(anglex1)

    angley1=angleyrad
    angley.append(angley1)

    Alist_in=np.zeros([numlayers,numfreqs])
    Alist_out=np.zeros(numlayers)

    for m in range(numlayers):
        anglextemp=np.zeros(numfreqs)
        angleytemp=np.zeros(numfreqs)
        nvectemp=np.zeros(numfreqs)
        avectemp=np.zeros(numfreqs)
        layertemp=Iso['layers'][m]
                    
        wout,aouttemp,nouttemp=layertemp.estimate(freqout)

        angleoutytemp=0.00
        angleoutxtemp=0.00
        koutz=0.00
        kouty=0.00
        koutx=0.00

        for i in range(numfreqs):
            anglex1temp,angley1temp=Angle(Iso,Las,m+1,i+1)
            w,a,n=layertemp.estimate(freqs[i])
            # NOTE: due to specific geometries used so far, it is unnecessary to have
            # r, theta, phi conversions
            #
            if (anglex1temp==0.00):
                anglez=np.pi/2-angley1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                kouty=2*np.pi*n*w*np.sin(angley1temp)*kcoeffs[i]+kouty
              
            else:
                anglez=np.pi/2-anglex1temp
                koutz=2*np.pi*n*w*np.sin(anglez)*kcoeffs[i]+koutz
                koutx=2*np.pi*n*w*np.sin(anglex1temp)*kcoeffs[i]+koutx

            anglextemp[i]=anglex1temp
            angleytemp[i]=angley1temp
            nvectemp[i]=n
            avectemp[i]=a

        angleoutxtemp=np.arctan(koutx/koutz)
        angleoutytemp=np.arctan(kouty/koutz)
        anglex.append(anglextemp)
        angley.append(angleytemp)      
        nvec.append(nvectemp)
        tk.append(layertemp['thickness'])
        angleouty.append(angleoutytemp)
        angleoutx.append(angleoutxtemp)
        nout.append(nouttemp)
        aout.append(aouttemp)

    launchanglex=np.arcsin(nouttemp*np.sin(angleoutxtemp))
    launchangley=np.arcsin(nouttemp*np.sin(angleoutytemp))

    angleoutx.append(launchanglex)
    angleouty.append(launchangley)

    for i in range(numfreqs):
        dttemp=float(0.00)
        for m in range(numlayers):
            anglextemp=anglex[m][i]
            angleytemp=angley[m][i]
            ntemp=nvec[m][i]
            atemp=avec[m][i]
            thick=tk[m]
            if (anglextemp==0):
                tkeff=thick/np.cos(angleytemp)
                abseff=atemp*tkeff
            else:
                tkeff=thick/np.cos(anglextemp)
                abseff=atemp*tkeff
            
            Alist_in[m][i]=abseff

    for m in range(numlayers):
        anglextemp=angleoutx[m][i]
        angleytemp=angleouty[m][i]
        ntemp=nout[m]
        thick=tk[m]
        if (anglextemp==0):
            tkeff=thick/np.cos(angleytemp)
            abseff=atemp*tkeff/np.log(10)
        else:
            tkeff=thick/np.cos(anglextemp)
            abseff=atemp*tkeff
        Alist_out[m]=abseff/np.log(10)

    return Alist_in, Alist_out


def applyabsorbances(Mlist, Alist_in, Alist_out=None):
    """Applies absorbances to the Mfactors calculated per layer in the IsoSample.

       Parameters
       ------
       Mlist:  1D array of M factors for the output FWM in each layer.     
       Alist_in: " 2D array with frequencies in 1st axis if (input) and layernums in second.  Elements are times (fsec)
       in which frequency i reaches the next layer with absorbances (log 10)
       Alist_out  1D array of output absorbances (log 10)
       Tdict:  dictionary of Fresnel coefficients with format described by that in Mcalc 
    
       Output
       -----
       Mlistnew:  Mlist scaled by Fresnel (if not None) and absorbance (if not None) losses
    
       """
    

    Mlistnew1=list()
    Mlistnew=list()

    if (Alist_out is None):
        Alist_out=list()
        Alistouttemp=float(0.000)
        for m in range(len(Alist_in)):
            Alist_out.append(float(Alistouttemp))

    numlayers=len(Alist_in)

    for m in range(numlayers):
        Mlistnewtemp=Mlist[m]
        for i in range(len(Alist_in[m])):
            Mlistnewtemp=10**(-Alist_in[m][i])*10**(-Alist_in[m][i])*Mlistnewtemp #note the squaring by the double 
            #multiplication as M factor is already a squared term
        Mlistnew1.append(Mlistnewtemp)   
    
    for m in range(numlayers):
        Mlistnewtemp=Mlistnew1[numlayers-1-m]
        Aouttemp=float(1.00)
        if ((numlayers-1-m) ==0):
            pass
        else:
            for n in range(numlayers-m+1,numlayers,1):
                Aouttemp=Aouttemp*10**(-Alist_out[n])  #this is NOT squared.
        Mlistnewtemp=Mlistnewtemp*Aouttemp
        Mlistnew.append(Mlistnewtemp)
    
    return Mlistnew


def applyfresneltrans(Mlist, Tdict=None):
    """Applies Fresnel coefficients to the Mfactors calculated per layer in the IsoSample.  Builds
        from all layers prior to that layer.  Transmission mode.

       Parameters
       ------
       Mlist_in:  2D array with frequencies in 1st axis if (input) and layernums in second.  Elements are times (fsec)
       in which frequency i  reaches the next layer.     
       Tdict:  dictionary of Fresnel coefficients with format described by that in Mcalc 
    
       Output
       -----
       Mlistnew:  Mlist scaled by Fresnel (if not None) losses
    
       """
    
    Mlistnew1=list()
    Mlistnew=list()

    if ((Tdict is None)):
        return Mlist

    Txin=Tdict['Txin']
    Tyin=Tdict['Tyin']
    Txout=Tdict['Txout']
    Tyout=Tdict['Tyout']
    numlayers=len(Mlist)

    for m in range(numlayers):
        Mlistnewtemp=Mlist[m]
        for i in range(len(Mlist[m])):
            for n in range(m): #builds up including previous layers
                Txintemp=Txin[n][i]
                Tyintemp=Tyin[n][i]
                if ((Txintemp==0.00) & (Tyintemp==0.00)):
                    Ttemp=0.00
                elif (Txintemp == 0.00):
                    Ttemp=Tyintemp
                elif (Txintemp == 0.00):
                    Ttemp=Txintemp
                else:
                    return ValueError("Unable to compute on input due to presence of x and y terms ")    
                Mlistnewtemp=Ttemp*Ttemp*Mlistnewtemp #again squared
        Mlistnew1.append(Mlistnewtemp)    
    
    for m in range(numlayers):
        Mlistnewtemp=Mlist[numlayers-m-1]
        Touttemp=float(1.00)
        for n in range(numlayers-m+1,numlayers,1):
            if (Txout[n]==0):
                Touttemp=Tyout[n]*Touttemp  #this is NOT squared.
            elif (Tyout[n]==0):
                Touttemp=Txout[n]*Touttemp
            else:
                return ValueError("Unable to compute on output due to presence of x and y terms")
        Mlistnewtemp=Mlistnewtemp*Touttemp
        Mlistnew.append(Mlistnewtemp)
    
    return Mlistnew
        

    