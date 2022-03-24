import numpy as np
from ._lasers import Lasers
from ._isosample import IsoSample
from ._isosample import Layer
from sympy import S, FiniteSet

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
    a tuple Mlist,Tdict consisting of:
    
    Mlist :  a complex array of phasemismatching factors for wavemixing at the output,
    currently only supporting four wave mixing models with supportedgeometries shown in the Laser object.
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

        for i in range(numfreqs):
            anglex1temp,angley1temp=Angle(Iso,Las,m+1,i+1)
            w,a,n=layertemp.estimate(freqs[i])
            koutz=0.00
            kouty=0.00
            koutx=0.00
            # NOTE: due to specific geometries used so far, it is unnecessary
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
                

#            angleoutxtemp=np.arctan(koutx/koutz)
#            angleoutytemp=np.arctan(kouty/koutz)

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
    launchanglexdeg=launchanglex*180/np.pi
    launchangley=np.arcsin(nouttemp*np.sin(angleoutytemp))
    launchangleydeg=launchangley*180/np.pi

    angleoutx.append(launchanglex)
    angleouty.append(launchangley)

    polout=_guessoutputpol(pols)
    
    Txin, Tyin, Txout, Tyout = _calculatetrans(nvec,anglex,angley,pols,nout,angleoutx,angleouty,polout)

    for m in range(numlayers):
        # Current code only utilizes projection along z.   In future other
        # projections can be incorporated.  (kx and ky are thus unused.)
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
        koutx=kout*np.sin(angleoutxtemp)
        kouty=kout*np.sin(angleoutytemp)
        koutz=np.sqrt(kout**2-koutx**2-kouty**2)

        dk=koutz-(kcoeffs[0]*kztemp[0]+kcoeffs[1]*kztemp[1]+kcoeffs[2]*kztemp[2])
        da=0.5*(aouttemp-(np.abs(kcoeffs[0])*avectemp[0]+np.abs(kcoeffs[1])*avectemp[1]+np.abs(kcoeffs[2])*avectemp[2]))
        Mc1=np.exp(-0.5*aouttemp*tktemp)
        Mc2=np.complex(np.cos(dk*tktemp),np.sin(dk*tktemp))*np.exp(-da*tktemp)
        Mc3=np.complex(-da*tktemp/((dk**2+da**2)*tktemp**2), -tktemp*dk/((dk**2+da**2)*tktemp**2))
        
        if (i==0):
            Mpre=1
        else:
            Mpre=Mc1*Mc2    

        Mctemp=Mpre*(Mc1)*(Mc2-1)*Mc3
        Mlist.append(Mctemp)

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
    return Mlist, Tdict


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
    (optional) frequency if not None replaces the frequency at that number with the given value
     prior to calculation

    Output
    ------
    tuple {anglex,angley} decomposed angles for the selected frequency in that layer according
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
    lasanglexin=lasanglex
    lasangleyin=lasangley

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


def EstimateAngle(Iso,Las,layernum,freqnum, frequency):
    '''
    Given an Isotropic Sample, a layer number and frequency number used in a geometry defined in Las object,
    determine the input angle (in air) required for phasematching that frequency.  Return NaN if
    a solution does not exist.
    
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
    Output
    ----
    Sympy: FiniteSet : theta = (deg) float of original angle in air needed for that PM condition.
           EmptySet if a solution cannot be found
           Reals if all real frequencies greater than zero are found.
    '''
    if (isinstance (Iso,IsoSample)== False):
        return ValueError("first argument not an object of class IsotropicSample")
    if (isinstance (Las,Lasers)== False):
        return ValueError("second argument not an object of class Lasers")
    if (freqnum < 1):
        return ValueError("freqnum cannot be less than 1")
    if (frequency < 0.00):
        return ValueError("frequency cannot be less than 0")

    freqs=Las.frequencies
    kcoeffs=Las.k_coeffs

    anglexrad=Las.anglesxrad
    angleyrad=Las.anglesyrad
    
    numfreqs=len(freqs)
    freqout=float(0.00)

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

    anglex1=anglexrad
    #anglex.append(anglex1)

    angley1=angleyrad
    #angley.append(angley1)
 
    freqs[freqnum-1]=frequency

    for i in range(numfreqs):
        freqout=freqout+kcoeffs[i]*freqs[i]

    flag=int(0)
    for m in range(layernum):
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

#                angleoutxtemp=np.arctan(koutx/koutz)
#                angleoutytemp=np.arctan(kouty/koutz)
            
            else:
                anglex1temp=0.00
                angley1temp=0.00
                anglez=0.00
            
            pass
            
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

        angleoutxtemp=np.arctan(koutx/koutz)
        angleoutytemp=np.arctan(kouty/koutz)

        if ((koutx==koutz) & (koutx==0.00)):
            if (m==(layernum-1)):
                flag=1
                break
        elif ((kouty==koutz) & (koutx==0.00)):
            if (m==(layernum-1)):
                flag=1
                break
        else:
            anglex.append(anglextemp)
            angley.append(angleytemp)      
            kx.append(kxtemp)
            ky.append(kytemp)
            kz.append(kztemp)
            avec.append(atemp)
            nvec.append(nvectemp)

            angleouty.append(angleoutytemp)
            angleoutx.append(angleoutxtemp)
            nout.append(nouttemp)

    if (flag==1):
        return S.Reals

    else:
        m = layernum-1
            # Current code only utilizes projection along z.   In future other
            # projections can be incorporated.  (kx and ky are thus unused.)

        kztemp=kz[m]
        kytemp=ky[m]
        kxtemp=kx[m]
    
        angleoutxtemp=angleoutx[m]
        angleoutytemp=angleouty[m]
        nouttemp=nout[m]
        nvectemp=nvec[m]

        kout=2*np.pi*freqout*nouttemp
        koutx=kout*np.sin(angleoutxtemp)
        kouty=kout*np.sin(angleoutytemp)
        koutz=np.sqrt(kout**2-koutx**2-kouty**2)

        dk=koutz-(kcoeffs[0]*kztemp[0]+kcoeffs[1]*kztemp[1]+kcoeffs[2]*kztemp[2]) # one of these is zero
        ksolve=dk
        theta=np.arcsin(ksolve/(2*np.pi*nvectemp[freqnum-1]*frequency*kcoeffs[freqnum-1]))
    
        if np.isnan(theta):
            #print("No real solution found, empty set returned")
            return S.EmptySet
        for m in range(layernum):
            n1=nvec[layernum-1-m][freqnum-1]
            if (layernum-1-m == 0):
                n2=1.0000
            else:
                n2=nvec[layernum-2-m][freqnum-1]
            theta=np.arcsin(n1/n2*np.sin(theta))
    
        thetadeg=theta/np.pi*180.00
        return FiniteSet(thetadeg)


def EstimateFrequency(Iso, Las, layernum, freqnum):
    '''Using the current frequency as first guess, solves for the nearest possible phasematching frequency
    at a fixed angle for that frequencynum in a given layer, using an iterative convergence.  Returns NaN
    if a solution cannot be found within the number of iterations internal to the procedure.
    
    Parameters
    ----
    Iso = The IsotropicSample object
    Las = The Lasers object with only supportedgeometry list capable of solutions
    layernum = layer in which to solve for angle
    freqnum = laser position (defined by geometry)

    Output
    ----
    Sympy: FiniteSet : frequency = (cm-1) float of frequency needed for that PM condition.
           EmptySet if a solution cannot be found.
           Reals if all real frequencies greater than zero are found.
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
    anglexrad=Las.anglesxrad
    angleyrad=Las.anglesyrad
    
    numfreqs=len(freqs)
    freqout=float(0.00)

    # These are 2D arrays where the 1st D is layer (1st and last are air) and 2nd are the input freqs
    anglex=list()
    angley=list()
    nvec=list()

    # These are 2D arrays where the 1st D is layer (no air layers) and 2nd are the input freqs
    kx=list()
    ky=list()
    kz=list()
      
    # These are 1D arrays where the D is layer (1st and last are air)
    nout=list()
    angleoutx=list()
    angleouty=list()

    anglex1=anglexrad
    #anglex.append(anglex1)

    angley1=angleyrad
    #angley.append(angley1)
 
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
                    
        wout,aouttemp,nouttemp=layertemp.estimate(freqout)
        angleoutytemp=0.00
        angleoutxtemp=0.00
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

#                angleoutxtemp=np.arctan(koutx/koutz)
#                angleoutytemp=np.arctan(kouty/koutz)
            
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

        angleoutxtemp=np.arctan(koutx/koutz)
        angleoutytemp=np.arctan(kouty/koutz)
        
        if ((koutx==koutz) & (koutx==0.00)):
            if (m==(layernum-1)):
                flag=1
                break
        elif ((kouty==koutz) & (koutx==0.00)):
            if (m==(layernum-1)):
                flag=1
                break
        else:       
            anglex.append(anglextemp)
            angley.append(angleytemp)      
            kx.append(kxtemp)
            ky.append(kytemp)
            kz.append(kztemp)

            nvec.append(nvectemp)

            angleouty.append(angleoutytemp)
            angleoutx.append(angleoutxtemp)
            nout.append(nouttemp)

    if (flag==1):
        return S.Reals
    else:
        m = layernum-1
        # Current code only utilizes projection along z.   In future other
        # projections can be incorporated.  (kx and ky are thus unused.)
        layertemp=Iso['layers'][m]

        kztemp=kz[m]
        kytemp=ky[m]
        kxtemp=kx[m]
        anglextemp=anglex[m]
        angleytemp=angley[m]

        angleoutxtemp=angleoutx[m]
        angleoutytemp=angleouty[m]
        nouttemp=nout[m]
        nvectemp=nvec[m]

        kout=2*np.pi*freqout*nouttemp
        koutx=kout*np.sin(angleoutxtemp)
        kouty=kout*np.sin(angleoutytemp)
        koutz=np.sqrt(kout**2-koutx**2-kouty**2)

        dk=koutz-(kcoeffs[0]*kztemp[0]+kcoeffs[1]*kztemp[1]+kcoeffs[2]*kztemp[2]) 
        maxiter=15
        tol=0.1
        wsolv1=freqs[freqnum-1]
        wsolv2=0
        i=0

        while(np.abs(wsolv1-wsolv2)>tol):
            wsolv2=wsolv1
            w,a,ntemp=layertemp.estimate(wsolv2)
            anglex1temp,angley1temp=Angle(Iso, Las, layernum, freqnum, frequency=wsolv1)
        
            if (anglex1temp != 0):
                tempangle=anglex1temp
            else:
                tempangle=angley1temp
            wsolv1=dk/(2*np.pi*ntemp*np.sin(tempangle)*kcoeffs[freqnum-1])
            i=i+1
            if (i >= maxiter):
                return S.EmptySet

        return FiniteSet(wsolv1)
