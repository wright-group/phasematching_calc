import numpy as np
from ._lasers import Lasers
from ._isosample import IsoSample
from ._isosample import Layer
from sympy import S, FiniteSet, Interval, oo
import time
import copy

# prototypes
Iso = IsoSample()
Las = Lasers()


def _array_compare(arr1, arr2):
    "compare two 2D arrays for equality, assumes each shape are same"
    booltest = True
    for n in range(len(arr1)):
        booltest = booltest & ((arr1[n] == arr2[n]).all())
    return booltest


def _calculate_trans(
    xmask, narr, anglexarr, angleyarr, polsvec, noutvec, angleoutxvec, angleoutyvec, polout
):
    """
    Uses Fresnel's equations to calculate transmission coefficients between layers in a sample,
    as well as after the last layer

    Assumes magnetic susceptibility of materials involved in layers to be zero.

    Parameters
    ---------
    xmask : list(int)
       "binary" masks actual float telling whether x or y operations are to be performed (ymask = NOT xmask)
    narr : list(float)
        2D array of n where the 1st D is layer and 2nd are the input freqs
    anglexarr : list(float)
        "  of angles in x (radians) where 1st D is first input angles, THEN number of layers
    angleyarr : list(float)
        "  of angles in y (radians) "
    noutvec : list(float)
        1D array of output n for all layers
    angleoutxvec : list(float)
        " of output angle in x for all layers AND additional output angle (assumed air)
    angleoutyvec : list(float)
        " of output angle in y for all layers "
    polsvec : list(int)
        array of polarizations for input lasers (int)...see list below for polout)
    polout : int
        output polarization (int) (based on tensor averaging, user must assume knowledge of
        isotropic tensors as calculation does not automatically return errors for an output polarization
        not expected to have signal due to isotropic averaging)
            1=vertical
            !1 = horizontal

    Return
    ------
    tuple:  Tin, Tout
    Tin= Transmission coefficient 2D array of 1D=frequencies, 2D=layers frequency in that layer according
    to the geometry and polarizations found in the Lasers object. Polarizations must be linear and array
    of the form 1 = Vertical and != 1 = Horizontal
    The normal list of geometries presume that each input
    is either x or y but not both.
    See calculate_absorbances for the calculation of absorbance in the layer itself.

    Tout= transmission coefficient through each layer for the output as well as the final exit into the
    air (m+1 layer array)
    """

    Tin = list()
    Tout = list()

    for m in range(len(noutvec)):
        if m == 0:
            n1vec = np.ones(len(polsvec))
        else:
            n1vec = narr[m - 1]
        n2vec = narr[m]
        anglex1vec = anglexarr[m]
        angley1vec = angleyarr[m]
        anglex2vec = anglexarr[m + 1]
        angley2vec = angleyarr[m + 1]

        Tintemp = np.zeros(len(polsvec))

        for i in range(len(polsvec)):
            n2 = n2vec[i]
            anglex2 = anglex2vec[i]
            angley2 = angley2vec[i]
            pol = polsvec[i]
            n1 = n1vec[i]
            anglex1 = anglex1vec[i]
            angley1 = angley1vec[i]

            if pol == 1:
                Rx = (n1 * np.cos(anglex1) - n2 * np.cos(anglex2)) ** 2 / (
                    n1 * np.cos(anglex1) + n2 * np.cos(anglex2)
                ) ** 2
                Ry = (n1 * np.cos(angley2) - n2 * np.cos(angley1)) ** 2 / (
                    n1 * np.cos(angley2) + n2 * np.cos(angley1)
                ) ** 2
                Tx = 1 - Rx
                Ty = 1 - Ry
            else:
                Rx = (n1 * np.cos(anglex2) - n2 * np.cos(anglex1)) ** 2 / (
                    n1 * np.cos(anglex2) + n2 * np.cos(anglex1)
                ) ** 2
                Ry = (n1 * np.cos(angley1) - n2 * np.cos(angley2)) ** 2 / (
                    n1 * np.cos(angley1) + n2 * np.cos(angley2)
                ) ** 2
                Tx = 1 - Rx
                Ty = 1 - Ry

            if xmask[i] == 0:
                Tintemp[i] = Ty
            else:
                Tintemp[i] = Tx

        Tin.append(Tintemp)

    i = len(xmask) - 1

    for m in range(len(noutvec)):
        n1 = noutvec[m]
        anglex2 = angleoutxvec[m + 1]
        angley2 = angleoutyvec[m + 1]
        anglex1 = angleoutxvec[m]
        angley1 = angleoutyvec[m]

        if m == (len(noutvec) - 1):
            n2 = 1.000
        else:
            n2 = noutvec[m + 1]

        if polout == 1:
            Rx = (n1 * np.cos(anglex1) - n2 * np.cos(anglex2)) ** 2 / (
                n1 * np.cos(anglex1) + n2 * np.cos(anglex2)
            ) ** 2
            Ry = (n1 * np.cos(angley2) - n2 * np.cos(angley1)) ** 2 / (
                n1 * np.cos(angley2) + n2 * np.cos(angley1)
            ) ** 2
            Tx = 1 - Rx
            Ty = 1 - Ry
        else:
            Rx = (n1 * np.cos(anglex2) - n2 * np.cos(anglex1)) ** 2 / (
                n1 * np.cos(anglex2) + n2 * np.cos(anglex1)
            ) ** 2
            Ry = (n1 * np.cos(angley1) - n2 * np.cos(angley2)) ** 2 / (
                n1 * np.cos(angley1) + n2 * np.cos(angley2)
            ) ** 2
            Tx = 1 - Rx
            Ty = 1 - Ry

        if xmask[i] == 0:
            Tout.append(Ty)
        else:
            Tout.append(Ty)

    return Tin, Tout


def _calculate_critical_angle(Iso, Las, layernum, freqnum, frequency=None):
    """Determines the critical angle for the layernum-1: layernum boundary.  Replaces value at freqnum
    by frequency if not None.

    Return
    -----
    tuple (n,crit) (float)
      refractive index, critical angle (radians) for the layernum-1:layernum boundary,
      or numpy.pi/2 if no critical angle found.
    """

    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")
    if freqnum < 1:
        return ValueError("freqnum cannot be less than 1")
    if layernum < 1:
        return ValueError("layernum cannot be less than 1")

    if frequency is not None:
        Las.frequencies[freqnum - 1] = frequency

    freq = Las.frequencies[freqnum - 1]

    m = layernum - 1
    if layernum == len(Iso["layers"]):
        w, a, nold = Iso["layers"][m].estimate(freq)
        n = float(1.000)
    else:
        w, a, nold = Iso["layers"][m].estimate(freq)
        w, a, n = Iso["layers"][m + 1].estimate(freq)

    if n / nold > 1:
        crit = np.pi / 2
    else:
        crit = np.arcsin(n / nold)

    return nold, crit


def _calculate_all_critical_angles(Iso, Las, freqnum=None, frequency=None):
    """Calculates all critical angles for all inputs at all layer boundaries.  Or returns
    numpy pi/2 if there is no critical angle at that combination.  Replaces freqnum with
    frequency if specified.

    Return
    -----
    out : list(float)
        2D arrray (1st D is layers m, 2nd is laser inputs i with float of the critical angle
        between m and m+1 for each element  or numpy.pi/2 if no critical angle found.
    """

    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")
    if freqnum < 1:
        return ValueError("freqnum cannot be less than 1")
    if frequency < 0.00:
        return ValueError("frequency cannot be less than 0")

    if (frequency is not None) & (freqnum is not None):
        Las.frequencies[freqnum - 1] = frequency

    mt = len(Iso["layers"])
    it = len(Las.frequencies)
    critarr = np.zeros([mt, it])
    narr = critarr

    for m in range(mt):
        for i in range(it):
            n, crit = _calculate_critical_angle(Iso, Las, m + 1, i + 1)
            critarr[m, i] = crit
            narr[m, i] = n

    return narr, critarr


def _calculate_internals(Iso, Las, zerofreq=False, zerofreqnum=1):
    """Returns a dictionary containing internal calculations used in many of the methods."""
    output = dict()
    freqs = Las.frequencies
    kcoeffs = Las.k_coeffs
    numlayers = len(Iso["layers"])
    anglexrad = Las.anglesxrad
    angleyrad = Las.anglesyrad
    xmask = Las.xmask

    numfreqs = len(freqs)
    freqout = float(0.00)
    for m in range(numfreqs):
        freqout = freqout + kcoeffs[m] * freqs[m]

    if freqout < 0:
        return ValueError("output frequency is less than zero")

    # These are 2D arrays where the 1st D is layer (1st and last are air) and 2nd are the input freqs
    anglex = list()
    angley = list()

    # These are 2D arrays where the 1st D is layer (no air layers) and 2nd are the input freqs
    kx = list()
    ky = list()
    kz = list()
    avec = list()
    nvec = list()
    # These are 1D arrays where the D is layer (1st and last are air)
    nout = list()
    angleoutx = list()
    angleouty = list()
    # These are 1D arrays where the D is layer (no air layers)
    tk = list()
    aout = list()

    anglex1 = anglexrad
    anglex.append(anglex1)

    angley1 = angleyrad
    angley.append(angley1)

    for m in range(numlayers):
        anglextemp = np.zeros(numfreqs)
        angleytemp = np.zeros(numfreqs)
        kxtemp = np.zeros(numfreqs)
        kytemp = np.zeros(numfreqs)
        kztemp = np.zeros(numfreqs)
        atemp = np.zeros(numfreqs)
        nvectemp = np.zeros(numfreqs)
        layertemp = Iso["layers"][m]

        wout, aouttemp, nouttemp = layertemp.estimate(freqout)

        angleoutytemp = 0.00
        angleoutxtemp = 0.00
        koutz = 0.00
        kouty = 0.00
        koutx = 0.00
        factor = 1.00
        for i in range(numfreqs):
            anglex1temp, angley1temp = angle(Iso, Las, m + 1, i + 1)
            w, a, n = layertemp.estimate(freqs[i])

            # NOTE: due to specific geometries used so far, it is unnecessary to have
            # r, theta, phi conversions
            #
            if xmask[i] == 0.00:
                anglez = np.pi / 2 - angley1temp
                koutz = 2 * np.pi * n * w * np.sin(anglez) * kcoeffs[i] + koutz
                kouty = 2 * np.pi * n * w * np.sin(angley1temp) * kcoeffs[i] + kouty
                anglex1temp = 0.00
            else:
                angley1temp = 0.00
                anglez = np.pi / 2 - anglex1temp
                koutz = 2 * np.pi * n * w * np.sin(anglez) * kcoeffs[i] + koutz
                koutx = 2 * np.pi * n * w * np.sin(anglex1temp) * kcoeffs[i] + koutx
            if zerofreq == True:
                if zerofreqnum is None:
                    return ValueError("zerofreqnum must be specified")
                else:
                    if i + 1 == zerofreqnum:
                        factor = 0.00
                    else:
                        factor = 1.00
            # this is put in as reminder that we are solving for this variable so the k's for this
            # one have to be set to zero if being used in a solver

            kx1 = 2 * np.pi * n * w * np.sin(anglex1temp) * factor
            ky1 = 2 * np.pi * n * w * np.sin(angley1temp) * factor
            kz1 = 2 * np.pi * n * w * np.sin(anglez) * factor
            anglextemp[i] = anglex1temp
            angleytemp[i] = angley1temp
            kxtemp[i] = kx1
            kytemp[i] = ky1
            kztemp[i] = kz1
            nvectemp[i] = n
            atemp[i] = a

        angleoutxtemp = np.arctan(koutx / koutz)
        angleoutytemp = np.arctan(kouty / koutz)
        anglex.append(anglextemp)
        angley.append(angleytemp)
        kx.append(kxtemp)
        ky.append(kytemp)
        kz.append(kztemp)
        avec.append(atemp)
        nvec.append(nvectemp)

        tk.append(layertemp["thickness"])
        aout.append(aouttemp)
        angleouty.append(angleoutytemp)
        angleoutx.append(angleoutxtemp)
        nout.append(nouttemp)

    output["anglex"] = anglex
    output["angley"] = angley
    output["angleoutx"] = angleoutx
    output["angleouty"] = angleouty
    output["n"] = nvec
    output["a"] = avec
    output["nout"] = nout
    output["aout"] = aout
    output["thicknesses"] = tk
    output["kx"] = kx
    output["ky"] = ky
    output["kz"] = kz
    output["freqout"] = freqout
    output["numlayers"] = numlayers
    output["numfreqs"] = numfreqs
    output["xmask"] = xmask
    output["kcoeffs"] = kcoeffs

    return output


def _guess_output_pol(polsvec):
    """Guesses output polarization based on isotropic averaging, given the input polarizations and that they are linear and
    aligned with either x or y but not a factor of both.  Returns a 1 or 0...1 means vertical, 0 means horizontal.
    """
    temppolsvec = np.zeros(len(polsvec))
    for i in range(len(polsvec)):
        if polsvec[i] != 1:
            temppolsvec[i] = 0
        else:
            temppolsvec[i] = 1
    s1 = np.sum(temppolsvec)
    if s1 == 3:
        return 1
    elif s1 == 2:
        return 0
    elif s1 == 1:
        return 1
    elif s1 == 0:
        return 0
    else:
        return ValueError(
            "Polarizations list not supported for current estimated output polarization"
        )


def _stack_critical_angles(narr, critarr, layernum, freqnum):
    """Up to a given layernum, finds the smallest of all critical angles of a freqnum and
    returns it.  This is necessary for the SolveAngle result when a large solution of angles are possible.

    Return
    -----
    crit : float
         angle in air (radians) to result in the smallest critangle in critarr up to layernum for freqnum.
    """

    critvec = list()
    nvec = list()

    for m in range(layernum):
        critvec.append(critarr[m][freqnum - 1])
        nvec.append(narr[m][freqnum - 1])

    crit = np.pi / 2

    for m in range(layernum):
        ncurrent = nvec[layernum - m - 1]
        critcurrent = critvec[layernum - m - 1]
        crit = np.minimum(crit, critcurrent)

        if m != layernum:
            nnew = nvec[layernum - m - 1]
        else:
            nnew = float(1.000)

        critnew = np.arcsin(nnew / ncurrent * np.sin(crit))

        if np.isnan(critnew):
            critnew = np.pi / 2

        crit = critnew

    return crit


def _m_plot(Iso, Las, layernum, freqnum, side=1):
    """Returns a full list of M factors for a complete input list of angles down to or up to zero.
    The parameter determines which of the two."""
    limit = int(75)  # degrees

    if side == 1:
        anglelist = list(range(limit, 0, -1))
    else:
        anglelist = list(range(-limit, 0, 1))

    mlist = list()
    Lastemp2 = copy.deepcopy(Las)
    Isotemp2 = copy.deepcopy(Iso)
    for k in range(layernum):
        Isotemp2.layers[k].suppress_absorbances()

    for m in range(len(anglelist)):
        Lastemp2.change_angle(freqnum, anglelist[m])
        Mfac, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
        # alist.append(m)
        mlist.append(Mfac[layernum - 1])

    return mlist


def calculate_original_crit_angle(Iso, Las, layernum, freqnum, frequency=None):
    """calculates the maximum original angle in air to result in the critical angle at the
    layernum : layernum+1 interface (or numpy pi/2 if there are no critical angles, or
    the angle in air of the smallest of any critical angles for that freqnum in any previous layers,
    converted to its original angle in air)

    Return
    ----
    angle: float
        (degrees) in air of critical angle
    """

    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")
    if freqnum < 1:
        return ValueError("freqnum cannot be less than 1")
    if layernum < 1:
        return ValueError("layernum cannot be less than 1")

    narr, critarr = _calculate_all_critical_angles(Iso, Las, freqnum, frequency)
    angle = _stack_critical_angles(narr, critarr, layernum, freqnum)
    angledeg = angle / np.pi * 180.00

    return angledeg


def launchangle(Iso, Las):
    """
    Calculates the launch angle in degrees for wavemixing for an IsotropicSample,

    Parameters
    ----
    Iso: an IsotropicSample object
    Las: a Lasers object

    Return
    ---
    tuple (angleoutx, angleouty) consisting of:

    angleoutx x-component of launchangle after the sample of the FWM output.
    angleouty y-component of launchangle after the sample of the FWM output.
    """
    try:
        output = _calculate_internals(Iso, Las, zerofreq=False, zerofreqnum=1)

        angleoutx = output["angleoutx"]
        angleouty = output["angleouty"]

        nout = output["nout"]

        numlayers = output["numlayers"]

        noutf = nout[numlayers - 1]
        angleoutxf = angleoutx[numlayers - 1]
        angleoutyf = angleouty[numlayers - 1]

        launchanglex = np.arcsin(noutf * np.sin(angleoutxf)) * 180 / np.pi
        launchangley = np.arcsin(noutf * np.sin(angleoutyf)) * 180 / np.pi

        return launchanglex, launchangley
    except:
        print("cannot calculate Launch Angle because output frequency is less than zero")
        return ValueError("output frequency is less than zero")


def m_calc(Iso, Las):
    """
    Calculates the phase mismatching factors for four wave mixing for an IsotropicSample, which can then
    be incorporated into plots or simulations (WrightSim)

    Parameters
    ----
    Iso: an IsotropicSample object
    Las: a Lasers object

    Return
    ---
    tuple (Mlist, Mlistdelta, tklist, Tdict) consisting of:

    Mlist :  list(float)
        a real array of phasemismatching factors for wavemixing at the output
        currently only supporting four wave mixing models with supportedgeometries shown in the Laser object.
    Mlistdelta : list (float)
        the phase (radians) of this term relative to terms from layers preceeding it, which can be used
        with the np.sqrt of Mlist to provide the Complex expression for the term inside the M^2 factor
        and can be used for small successive layer calculations
    tklist :  list(float)
        the effective thickness of each layer as pertaining to the launched output wave, i.e.
        layer thickness / cosine(angleout)
    Tdict :  dict
        dictionary with entries related to transmission coefficients of the lasers and output
        through the sample layers.  NOTE:  This application has not been performed on the entries in Mlist.
        The coefficients rely on linear polarizations defined in the Laser object.  User must perform
        the multiplication.  Methods to do so are available in this module.  See "applyfresneltrans".

        Keys in the Tdict are:
            Tdict['Tin'] : list(float)  (2Darr) Fresnel transmission coefficients of the laser inputs based on the geometry,
                with 1st D as layer and 2nd as freqnum

            Tdict['Tout'] : list(float)  (1D) Fresnel transmission coefficients of  the output based on the geometry .

            Tdict['launchangledeg'] : list(float) 1D launch angle of output in degrees in air at cartesian coordinate
                specified by geometry
    """
    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")

    pols = Las.polarizations

    Mlist = list()
    tklist = list()
    Mphaselist = list()

    try:
        output = _calculate_internals(Iso, Las, zerofreq=False, zerofreqnum=1)

        anglex = output["anglex"]
        angley = output["angley"]
        angleoutx = output["angleoutx"]
        angleouty = output["angleouty"]
        nvec = output["n"]
        avec = output["a"]
        nout = output["nout"]
        aout = output["aout"]
        tk = output["thicknesses"]
        kx = output["kx"]
        ky = output["ky"]
        kz = output["kz"]
        freqout = output["freqout"]
        numlayers = output["numlayers"]
        numfreqs = output["numfreqs"]
        xmask = output["xmask"]
        kcoeffs = output["kcoeffs"]

        noutf = nout[numlayers - 1]
        angleoutxf = angleoutx[numlayers - 1]
        angleoutyf = angleouty[numlayers - 1]

        launchanglex = np.arcsin(noutf * np.sin(angleoutxf))
        launchangley = np.arcsin(noutf * np.sin(angleoutyf))

        angleoutx.append(launchanglex)
        angleouty.append(launchangley)

        tkeffvec = np.zeros(numfreqs)

        polout = _guess_output_pol(pols)

        Tin, Tout = _calculate_trans(
            xmask, nvec, anglex, angley, pols, nout, angleoutx, angleouty, polout
        )

        Mctemp = 1

        koutx = 0.00
        kouty = 0.00
        koutz = 0.00

        n = len(xmask) - 1
        for m in range(numlayers):
            avectemp = avec[m]
            kztemp = kz[m]
            kytemp = ky[m]
            kxtemp = kx[m]

            tktemp = tk[m]
            anglextemp = anglex[m + 1]  ## NOTE!! +1
            angleytemp = angley[m + 1]  ##
            aouttemp = aout[m]
            angleoutxtemp = angleoutx[m]
            angleoutytemp = angleouty[m]
            nouttemp = nout[m]

            if xmask[n] == 0.00:
                tkeff = tktemp / np.cos(angleoutytemp)
                kout = 2.00 * np.pi * freqout * nouttemp * tkeff
                koutx = 0.00
                koutz = kout * np.cos(angleoutytemp)
                kouty = kout * np.sin(angleoutytemp)
                ko = 2.00 * np.pi * freqout * nouttemp
                kox = 0.00
                koz = ko * np.cos(angleoutytemp)
                koy = ko * np.sin(angleoutytemp)
            else:
                tkeff = tktemp / np.cos(angleoutxtemp)
                kout = 2.00 * np.pi * freqout * nouttemp * tkeff
                kouty = 0.00
                koutz = kout * np.cos(angleoutxtemp)
                koutx = kout * np.sin(angleoutxtemp)
                ko = 2.00 * np.pi * freqout * nouttemp
                koy = 0.00
                koz = ko * np.cos(angleoutxtemp)
                kox = ko * np.sin(angleoutxtemp)

            ksumx = ksumy = ksumz = 0.0000
            # kz=0.000

            for i in range(numfreqs):
                if xmask[i] == 0.00:
                    tkeffvec[i] = tktemp / np.cos(angleytemp[i])
                    k = np.sqrt(kytemp[i] ** 2 + kztemp[i] ** 2)
                    kproj = k * np.cos(angleytemp[i] - angleoutytemp)
                else:
                    tkeffvec[i] = tktemp / np.cos(anglextemp[i])
                    k = np.sqrt(kxtemp[i] ** 2 + kztemp[i] ** 2)
                    kproj = k * np.cos(anglextemp[i] - angleoutxtemp)
                """
                #ksumx = (kcoeffs[i] * kxtemp[i]) * tkeffvec[i] + ksumx
                #ksumy = (kcoeffs[i] * kytemp[i]) * tkeffvec[i] + ksumy
                #ksumz = (kcoeffs[i] * kztemp[i]) * tkeffvec[i] + ksumz
                #ksumx = (kcoeffs[i] * kxtemp[i]) * tkeff + ksumx #B
                #ksumy = (kcoeffs[i] * kytemp[i]) * tkeff + ksumy #B
                #ksumz = (kcoeffs[i] * kztemp[i]) * tkeff + ksumz #B
                """
                ksumx = 0  # C
                ksumy = 0
                ksumz = (kcoeffs[i] * kproj) * tkeff + ksumz
                # kz = (kcoeffs[i] * kproj) + ksumz

            k4 = np.sqrt(ksumx**2 + ksumy**2 + ksumz**2)
            dkl = k4 - kout

            # dk = kz-ko

            dal = aouttemp * tkeff * 0.5
            da = aouttemp * 0.5

            for i in range(numfreqs):
                dal = dal - (np.abs(kcoeffs[i]) * avectemp[i] * tkeffvec[i]) * 0.5  # C, #A
                da = da - (np.abs(kcoeffs[i]) * avectemp[i]) * 0.5  # C, #A
                # dal = dal - 0.5 * (np.abs(kcoeffs[i]) * avectemp[i] * tkeff)

            Mc1 = np.exp(-aouttemp * tkeff)

            if (dal == 0.00) & (dkl == 0.00):  # limiting case
                Mc2 = 1.00
                Mphasedelta = 0.00
                # Mphasedelta=np.pi/2
            else:

                Mc2 = ((1 - np.exp(dal)) ** 2 + 4 * np.exp(dal) * (np.sin(dkl / 2)) ** 2) / (
                    dal**2 + (dkl) ** 2
                )
                """
                Mc2 = (1+ np.exp(2*dal)-2*np.exp(dal)*np.cos(dkl)) / (
                    dal**2 + (dkl) ** 2
                )
                """

                Mphasedelta = np.arctan(
                    (dkl + np.exp(dal) * (-dkl * np.cos(dkl) + dal * np.sin(dkl)))
                    / (-dal + np.exp(dal) * (dal * np.cos(dkl) + dkl * np.sin(dkl)))
                )

                """
                Mphasedelta = np.arctan(
                    (-dal + np.exp(dal) * (dkl * np.sin(dkl) + dal * np.cos(dkl)))
                    / (dkl + np.exp(dal) * (-dkl * np.cos(dkl) + dal * np.sin(dkl)))
                ) """  # old

            Mctemp = Mc1 * Mc2
            Mlist.append(Mctemp)
            tklist.append(tkeff)
            Mphaselist.append(Mphasedelta)
        # angleoutxtemp and angleoutytemp can be converted to launch angles via Snell's Law to calculate
        # launch angle in air afterwards
        if xmask[n] == 0:
            launchangledeg = np.arcsin(nouttemp * np.sin(angleoutytemp)) * 180.00 / np.pi
        else:
            launchangledeg = np.arcsin(nouttemp * np.sin(angleoutxtemp)) * 180.00 / np.pi

        Tdict = dict()
        Tdict["Tout"] = Tout
        Tdict["Tin"] = Tin
        Tdict["launchangledeg"] = launchangledeg
        return Mlist, Mphaselist, tklist, Tdict
    except:
        print("cannot calculate M factor because output frequency is less than zero")
        return ValueError("output frequency is less than zero")


def angle(Iso, Las, layernum, freqnum, frequency=None):
    """
    Uses Snell's Law to calculate the angle a laser makes in a specific layer of a sample,
    the laser defined by freqnum and layer by layernum

    Parameters
    ---------
    Iso : IsoSample object
    Las : Lasers object
    layernum : int
       layer number
    freqnum : int
       frequency number (index number)
    frequency : float (optional)
       if not None replaces the frequency at that number with the given value prior to calculation

    Return
    ------
    tuple: anglex,angley
        anglex, angley: float
        radians) for the selected frequency in that layer according
        to the geometry found in the Lasers object.  NOTE:  User must determine which (or both)
        are to be used.  The xmask or ymask of a geometry is used in methods in this module to do so.
    """
    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")
    if freqnum < 1:
        return ValueError("freqnum cannot be less than 1")
    if layernum < 1:
        return ValueError("layernum cannot be less than 1")

    lasfreq = Las.frequencies[freqnum - 1]
    lasanglex = Las.anglesxrad[freqnum - 1]
    lasangley = Las.anglesyrad[freqnum - 1]

    if frequency is not None:
        lasfreq = frequency

    for m in range(layernum):
        if m == 0:
            nold = float(1.000)
            aold = float(0.000)
        else:
            nold = n
            aold = a
        w, a, n = Iso["layers"][m].estimate(lasfreq)
        lasanglex = np.arcsin(nold / n * np.sin(lasanglex))
        lasangley = np.arcsin(nold / n * np.sin(lasangley))
    return lasanglex, lasangley

    """
    w, a, n = Iso["layers"][layernum-1].estimate(lasfreq)
    lasanglex = np.arcsin(1.00/ n * np.sin(lasanglex))
    lasangley = np.arcsin(1.00 / n * np.sin(lasangley))
    return lasanglex, lasangley
    """


def solve_angle(Iso, Las, layernum, freqnum, frequency=None, isclose=False, amt=None):
    """
    Given an Isotropic Sample, a layer number and frequency number used in a geometry defined in Las object,
    determine the input angle (in air) required for phasematching that frequency. Uses Sympy Set.
    Returns an empty FiniteSet if a solution cannot be found within an internal convergence iteration series,
    usually meaning there is no solution.

    Parameters
    ----
    Iso : The IsotropicSample object
    Las : The Lasers object
         with only supportedgeometry list capable of solutions
    layernum : int
         layer in which to solve for angle
    freqnum : int
         laser position to solve for (defined by geometry)
    frequency : float (optional)
         frequency at that laser position rather than the value indicated in the Lasers object, if desired.
    isclose : bool (optional)
         if set to True it will only search for a solution close to the angle shown in the Lasers object
            with a faster algorithm

    All other frequencies, k coefficients, refractive indexes are taken from the objects.

    Return
    ----
    tuple:  Sympy, Amount

        Sympy: FiniteSet : {theta1,theta2}
           theta1,theta2: float
              original angles in air (degrees) needed for that PM condition.
              theta1 and/or theta2 may be missing if either or both cannot be found as solution
              if isclose is set to True, only one solution at most will be found.
           OR Interval(0,anglemax)
              anglemax (deg) if all angles allowed until restricted by a critical angle before the layernum

        Amount: float
            passes out the most recent amount of changed used per step in the solver.  Useful for
            speeding up the calculation of an "isclose" routine on repeated uses.

    """
    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")
    if freqnum < 1:
        return ValueError("freqnum cannot be less than 1")

    freqs = Las.frequencies
    Lastemp = copy.deepcopy(Las)
    Isotemp = copy.deepcopy(Iso)
    Lastemp2 = copy.deepcopy(Las)
    Isotemp2 = copy.deepcopy(Iso)

    numfreqs = len(freqs)

    flag = int(0)
    flag2 = int(1)
    flag3 = int(1)

    if frequency is not None:
        freqs[freqnum - 1] = frequency
        if frequency < 0.00:
            return ValueError("frequency cannot be less than 0")
        Lastemp.change_freq(freqnum, frequency)

    try:
        output = _calculate_internals(Isotemp, Lastemp, zerofreq=True, zerofreqnum=freqnum)

        kx = output["kx"]
        ky = output["ky"]
        kz = output["kz"]
        kcoeffs = output["kcoeffs"]
        kxtemp = kx[layernum - 1]
        kytemp = ky[layernum - 1]
        kztemp = kz[layernum - 1]
        ksumx = ksumy = ksumz = 0.000

        for i in range(numfreqs):
            ksumx = kcoeffs[i] * kxtemp[i] + ksumx
            ksumy = kcoeffs[i] * kytemp[i] + ksumy
            ksumz = kcoeffs[i] * kztemp[i] + ksumz
        if (ksumx == ksumz) & (ksumx == 0.00):
            flag = 1
        elif (ksumy == ksumz) & (ksumx == 0.00):
            flag = 1

        if flag == 1:
            angledeg = calculate_original_crit_angle(Iso, Las, layernum, freqnum, frequency)
            return Interval(0, angledeg), 0.00
        else:
            m = layernum - 1
            tol = 0.001
            iter = 80
            for k in range(layernum):
                Isotemp2.layers[k].suppress_absorbances()

            if isclose:
                if amt is None:
                    amt = 0.05

                Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
                angle = Lastemp2.anglesairdeg[freqnum - 1]
                magMtestc = np.abs(Mtest[m])

                anglep = Lastemp2.anglesairdeg[freqnum - 1] + amt * 1.00
                Lastemp.change_angle(freqnum, anglep)
                Mtestp, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp)
                magMtestp = np.abs(Mtestp[m])

                anglen = Lastemp2.anglesairdeg[freqnum - 1] + 2 * amt * -1.00
                Lastemp.change_angle(freqnum, anglen)
                Mtestn, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp)
                magMtestn = np.abs(Mtestn[m])

                errorp = magMtestp - magMtestc
                errorn = magMtestn - magMtestc

                if errorp > errorn:
                    dir = 1.00
                    magMtest = magMtestc
                else:
                    dir = -1.00
                    magMtest = magMtestc

                error2 = 1.000 - magMtestc
                error1 = error2
                b = 0
                while error2 > tol:
                    b = b + 1
                    if error2 > error1:
                        dir = (-1.00) * dir
                    if np.abs(error2) < 0.33 * np.abs(error1):
                        amt = 0.1 * amt
                    error1 = error2
                    angle = Lastemp2.anglesairdeg[freqnum - 1] + amt * dir
                    Lastemp2.change_angle(freqnum, angle)
                    Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
                    magMtest = np.abs(Mtest[m])
                    error2 = 1 - magMtest
                    if b > iter:
                        flag2 = 1
                        break
                if np.isclose(magMtest, 1.00, rtol=tol * 10):
                    flag2 = 0
                angle2 = float("nan")
                flag3 = 1
            else:  # if isclose is false, both (if possible) solutions are located (two while loops)
                mlist = _m_plot(Isotemp2, Lastemp2, layernum, freqnum, side=1)
                mlist2 = _m_plot(Isotemp2, Lastemp2, layernum, freqnum, side=-1)
                mlist.reverse()
                max1 = max(mlist)
                max1ind = mlist.index(max(mlist))
                max2 = max(mlist2)
                max2ind = mlist2.index(max(mlist2)) - 75

                dir = 1.000
                amt = 0.5
                Lastemp2.change_angle(freqnum, max1ind)
                Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
                angle = Lastemp2.anglesairdeg[freqnum - 1]
                magMtest1 = np.abs(Mtest[m])
                error2 = 1.000 - magMtest1
                error1 = error2
                b = 0
                while error2 > tol:
                    b = b + 1
                    if error2 > error1:
                        dir = (-1.00) * dir
                    if np.abs(error2) < 0.33 * np.abs(error1):
                        amt = 0.1 * amt
                    error1 = error2
                    angle = Lastemp2.anglesairdeg[freqnum - 1] + amt * dir
                    Lastemp2.change_angle(freqnum, angle)
                    Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
                    magMtest1 = np.abs(Mtest[m])
                    error2 = 1 - magMtest1
                    if b > iter:
                        flag2 = 1
                        break

                dir = 1.00
                amt = 0.5
                Lastemp2.change_angle(freqnum, max2ind)
                Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
                angle2 = Lastemp2.anglesairdeg[freqnum - 1]
                magMtest2 = np.abs(Mtest[m])
                error2 = 1.000 - magMtest2
                error1 = error2
                b = 0
                while error2 > tol:
                    b = b + 1
                    if error2 > error1:
                        dir = (-1.00) * dir
                    if np.abs(error2) < 0.33 * np.abs(error1):
                        amt = 0.1 * amt
                    error1 = error2
                    angle2 = Lastemp2.anglesairdeg[freqnum - 1] + amt * dir
                    Lastemp2.change_angle(freqnum, angle2)
                    Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp2)
                    magMtest2 = np.abs(Mtest[m])
                    error2 = 1 - magMtest2
                    if b > iter:
                        flag3 = 1
                        break

                if np.isclose(magMtest1, 1.00, rtol=tol * 10):
                    flag2 = 0
                if np.isclose(magMtest2, 1.00, rtol=tol * 10):
                    flag3 = 0

            if (flag2 == 1) & (flag3 == 1):
                return FiniteSet(), amt
            elif flag2 == 1:
                return FiniteSet(angle2), amt
            elif flag3 == 1:
                return FiniteSet(angle), amt
            else:
                return FiniteSet(angle, angle2), amt
    except:
        print("cannot solve angle because output frequency is less than zero")
        return ValueError("output frequency is less than zero")


def solve_frequency(Iso, Las, layernum, freqnum, amt=None, isclose=False):
    """Using the current frequency as first guess, solves for the nearest possible phasematching frequency
    at a fixed angle for that frequencynum in a given layer, using an iterative convergence.  Uses Sympy Set.
    Returns an empty FiniteSet if a solution cannot be found within an internal convergence iteration series,
    usually meaning there is no solution.

    Parameters
    ----
    Iso : The IsotropicSample object
    Las : The Lasers object
       with only supportedgeometry list capable of solutions
    layernum : int
       layer number in which to solve for angle
    freqnum : int
       laser position to solve for(defined by geometry)
    amt : float (optional)
       amount to change frequency by per convergence step.  Approximates internally if set to None.
    isclose : bool (optional)
       uses finer algorithm if set to True

    Return
    ----
    tuple:  Sympy, amount
        Sympy: FiniteSet : {frequency} : frequency float
             (cm-1) of frequency needed for that PM condition.
             OR FiniteSet {} if a solution cannot be found.
             OR Interval(0,oo) if all real frequencies greater than zero are found
        amount:  float
            Last amount of frequency per convergence step.  Useful as a return for future solves (see Inputs)

    Interval is not restricted to critical angle assessments at this time.
    """
    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")
    if freqnum < 1:
        return ValueError("freqnum cannot be less than 1")
    if layernum < 1:
        return ValueError("layernum cannot be less than 1")

    freqs = Las.frequencies
    flag = int(0)
    numfreqs = len(freqs)

    Isotemp = copy.deepcopy(Iso)
    Isotemp2 = copy.deepcopy(Iso)
    Lastemp = copy.deepcopy(Las)

    if isclose:
        if amt is None:
            amt = 0.001 * freqs[freqnum - 1]
            tol = 0.0001
        else:

            tol = 0.0001
    else:
        if amt is None:
            amt = 0.01 * freqs[freqnum - 1]  # a guess
            tol = 0.001
        else:

            tol = 0.001

    try:
        output = _calculate_internals(Isotemp, Lastemp, zerofreq=True, zerofreqnum=freqnum)

        kx = output["kx"]
        ky = output["ky"]
        kz = output["kz"]
        numfreqs = output["numfreqs"]
        kcoeffs = output["kcoeffs"]

        ksumx = ksumy = ksumz = 0.000

        kxtemp = kx[layernum - 1]
        kytemp = ky[layernum - 1]
        kztemp = kz[layernum - 1]

        for i in range(numfreqs):
            ksumx = kcoeffs[i] * kxtemp[i] + ksumx
            ksumy = kcoeffs[i] * kytemp[i] + ksumy
            ksumz = kcoeffs[i] * kztemp[i] + ksumz

        if (ksumx == ksumz) & (ksumx == 0.00):
            flag = 1
        elif (ksumy == ksumz) & (ksumx == 0.00):
            flag = 1

        if flag == 1:
            return Interval(0, oo)  # currently does not attempt a walkback of high frequencies
            # at the given angle to see if it is reflected at a critical angle
        else:
            m = layernum - 1

            for k in range(layernum):
                Isotemp2.layers[k].suppress_absorbances()

            dir = 1.00
            Mtestc, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp)
            magMtestc = np.abs(Mtestc[m])

            freq = Lastemp.frequencies[freqnum - 1] + amt * dir
            Lastemp.change_freq(freqnum, freq)
            Mtestp, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp)
            magMtestp = np.abs(Mtestp[m])

            freq = Lastemp.frequencies[freqnum - 1] - 2 * amt * dir
            Lastemp.change_freq(freqnum, freq)
            Mtestn, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp)
            magMtestn = np.abs(Mtestn[m])

            errorp = magMtestp - magMtestc
            errorn = magMtestn - magMtestc

            if errorp > errorn:
                dir = 1.00
                magMtest = magMtestc
            else:
                dir = -1.00
                magMtest = magMtestc

            error1 = 1 - magMtest
            error2 = error1
            iter = 500
            b = 0
            while error2 > tol:
                b = b + 1
                if error2 > error1:
                    dir = (-1.00) * dir
                if np.abs(error2) < 0.33 * np.abs(error1):  # the 0.33 is a guess
                    amt = 0.1 * amt
                error1 = error2
                freq = Lastemp.frequencies[freqnum - 1] + amt * dir
                Lastemp.change_freq(freqnum, freq)
                Mtest, Mdelta, tklist, Tdict = m_calc(Isotemp2, Lastemp)
                magMtest = np.abs(Mtest[m])
                error2 = 1 - magMtest
                if b > iter:
                    flag = 2
                    break

            if flag == 2:
                return FiniteSet(), amt
            else:
                return FiniteSet(freq), amt
    except:
        print("cannot solve frequency because output frequency is less than zero")
        return ValueError("output frequency is less than zero")


def calculate_ts(Iso, Las):
    """Calculate the times each pulse makes with respect to the first pulse (first on the Las.frequencies
       list)  as each passes through the IsoSample.

       Uses the refractive indexes, angles of incidence, and thickness per sample, to determine the change in
       timing each pulse makes through a sample relative to the first pulse assuming all are perfectly overlapped
       before entering the IsoSample.

       Is used to confirm that the pulse overlap in thicker samples is satisfactory or not satisfactory enough
       to determine if delay changes are responsible for a limit to the interaction in such a sample.

    Parameters
    -----
    Iso :  The IsoSample object
    Las : The Lasers object

    Return
    -----
    tuple : (t_chart_in, t_chart_out)
            t_chart_in: list(float)
                2D array with frequencies in 1st axis and layernums in second.  Elements are times (fsec)
                in which frequency i  reaches the next layer.
            t_chart_out: list(float)
                1D array with layernums for the output at the given kcoeffs.
    """

    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")

    numlayers = len(Iso["layers"])
    freqs = Las.frequencies
    cvac = 29979245800  # cm/sec
    numfreqs = len(freqs)
    t_chart_in = np.zeros([numlayers, numfreqs])
    t_chart_out = np.zeros(numlayers)

    try:
        output = _calculate_internals(Iso, Las, zerofreq=False, zerofreqnum=1)

        anglex = output["anglex"]
        angley = output["angley"]
        angleoutx = output["angleoutx"]
        angleouty = output["angleouty"]
        nvec = output["n"]
        nout = output["nout"]
        tk = output["thicknesses"]
        numlayers = output["numlayers"]
        numfreqs = output["numfreqs"]
        xmask = output["xmask"]

        noutf = nout[numlayers - 1]
        angleoutxf = angleoutx[numlayers - 1]
        angleoutyf = angleouty[numlayers - 1]

        launchanglex = np.arcsin(noutf * np.sin(angleoutxf))
        launchangley = np.arcsin(noutf * np.sin(angleoutyf))

        angleoutx.append(launchanglex)
        angleouty.append(launchangley)

        for i in range(numfreqs):
            dttemp = float(0.00)
            for m in range(numlayers):
                anglextemp = anglex[m + 1][i]  ##NOTE the +1
                angleytemp = angley[m + 1][i]  ##
                ntemp = nvec[m][i]
                thick = tk[m]
                if xmask[i] == 0:
                    tkeff = thick / np.cos(angleytemp)
                else:
                    tkeff = thick / np.cos(anglextemp)
                dttemp = dttemp + tkeff * ntemp / cvac * 1e15
                t_chart_in[m][i] = dttemp

        dttemp = float(0.00)
        for m in range(numlayers):
            anglextemp = angleoutx[m]
            angleytemp = angleouty[m]
            ntemp = nout[m]
            thick = tk[m]
            if xmask[numfreqs] == 0:
                tkeff = thick / np.cos(angleytemp)
            else:
                tkeff = thick / np.cos(anglextemp)
            dttemp = dttemp + tkeff * ntemp / cvac * 1e15
            t_chart_out[m] = dttemp

        return t_chart_in, t_chart_out
    except:
        print("cannot calculate delta ts because output frequency is less than zero")
        return ValueError("output frequency is less than zero")


def calculate_absorbances(Iso, Las):
    """Calculate the absorbances each laser (including output) makes per layer in the IsoSample.

     Parameters
     -----
     Iso:  The IsoSample object
     Las:  The Lasers object

    Return
    ------
    tuple(Alist_in,Alist_out):
         Alist_in: list(float)
             2D array with frequencies in 1st axis and layernums in second.  Elements are log10 absorbances
             in which frequency i reaches the next layer.
         Alist_out: list(float)
             1D array of log10 absorbances at the output frequency as it reaches the next layer.
    """

    if isinstance(Iso, IsoSample) == False:
        return ValueError("first argument not an object of class IsotropicSample")
    if isinstance(Las, Lasers) == False:
        return ValueError("second argument not an object of class Lasers")

    try:
        output = _calculate_internals(Iso, Las, zerofreq=False, zerofreqnum=1)
    except:
        print("cannot calculate absorbances because output frequency is less than zero")
        return ValueError("output frequency is less than zero")

    anglex = output["anglex"]
    angley = output["angley"]
    angleoutx = output["angleoutx"]
    angleouty = output["angleouty"]
    avec = output["a"]
    nout = output["nout"]
    tk = output["thicknesses"]
    numlayers = output["numlayers"]
    numfreqs = output["numfreqs"]
    xmask = output["xmask"]

    noutf = nout[numlayers - 1]
    angleoutxf = angleoutx[numlayers - 1]
    angleoutyf = angleouty[numlayers - 1]

    launchanglex = np.arcsin(noutf * np.sin(angleoutxf))
    launchangley = np.arcsin(noutf * np.sin(angleoutyf))

    angleoutx.append(launchanglex)
    angleouty.append(launchangley)

    Alist_in = np.zeros([numlayers, numfreqs])
    Alist_out = np.zeros(numlayers)
    n = len(xmask) - 1

    for i in range(numfreqs):
        for m in range(numlayers):
            anglextemp = anglex[m + 1][i]  ##NOTE the +1
            angleytemp = angley[m + 1][i]  ##
            atemp = avec[m][i]
            thick = tk[m]
            if xmask[i] == 0:
                tkeff = thick / np.cos(angleytemp)
                abseff = atemp * tkeff
            else:
                tkeff = thick / np.cos(anglextemp)
                abseff = atemp * tkeff
            Alist_in[m][i] = abseff

    for m in range(numlayers):
        anglextemp = angleoutx[m]
        angleytemp = angleouty[m]
        thick = tk[m]
        if xmask[n] == 0:
            tkeff = thick / np.cos(angleytemp)
            abseff = atemp * tkeff / np.log(10)
        else:
            tkeff = thick / np.cos(anglextemp)
            abseff = atemp * tkeff
        Alist_out[m] = abseff / np.log(10)

    return Alist_in, Alist_out


def apply_absorbances(Mlist, Alist_in, Alist_out=None):
    """Applies absorbances to the Mfactors calculated per layer in the IsoSample, based on those absorbances preceding the layer.

    Parameters
    ------
    Mlist :  list(float)
         1D array of M factors for the output FWM in each layer.
    Alist_in : list(float)
         " 2D array with frequencies in 1st axis if (input) and layernums in second.  Elements are log10 absorbances
         in which frequency i reaches the next layer
    Alist_out : list(float)
         1D array of output absorbances (log 10)
    Tdict:  dictionary of Fresnel coefficients with format described by that in Mcalc

    Return
    -----
    Mlistnew:  list(float)
         Mlist scaled by Fresnel (if not None) and absorbance (if not None) losses

    See m_calc for description of Mlist and calculate_absorbances for more on the Alists.  See m_calc for
    description of Tdict.

    """

    Mlistnew1 = list()
    Mlistnew = list()

    if Alist_out is None:
        Alist_out = list()
        Alistouttemp = float(0.000)
        for m in range(len(Alist_in)):
            Alist_out.append(float(Alistouttemp))

    numlayers = len(Alist_in)

    for m in range(numlayers):
        Mlistnewtemp = Mlist[m]
        for i in range(len(Alist_in[m])):
            for n in range(m):
                if m == 0:  # technically not needed as range(0) leads to a null but kept in
                    pass
                else:
                    Mlistnewtemp = (
                        10 ** (-Alist_in[n][i]) * 10 ** (-Alist_in[n][i]) * Mlistnewtemp
                    )  # note the squaring by the double
            # multiplication as M factor is already a squared term
        Mlistnew1.append(Mlistnewtemp)

    for m in range(numlayers):
        Mlistnewtemp = Mlistnew1[numlayers - 1 - m]
        Aouttemp = float(0.00)
        if (numlayers - 1 - m) == 0:
            pass
        else:
            for n in range(numlayers - 1, numlayers - m - 1, -1):
                if m == 0:
                    pass
                else:
                    Aouttemp = Aouttemp + Alist_out[n]  # this is NOT squared.
        Mlistnewtemp = Mlistnewtemp * 10 ** (-Aouttemp)
        Mlistnew.append(Mlistnewtemp)

    Mlistnew.reverse()
    return Mlistnew


def apply_trans(Mlist, Tdict=None):
    """Applies Fresnel coefficients to the Mfactors calculated per layer in the IsoSample.  Builds
     from all layers prior to that layer.  This is for a strictly transmission mode experiment.

    Parameters
    ------
    Mlist_in : list(float)
         2D array with frequencies in 1st axis if (input) and layernums in second.  Elements are times (fsec)
         in which frequency i  reaches the next layer.
    Tdict:  dict
         dictionary of Fresnel coefficients with format described by that in m_calc

    Return
    -----
    Mlistnew: list(float)
         Mlist scaled by Fresnel (if not None) losses

     See m_calc for description of Mlist and calculate_absorbances for more on the Alists.  See m_calc for
    description of Tdict.
    """

    Mlistnew1 = list()
    Mlistnew = list()

    if Tdict is None:
        return Mlist

    Tin = Tdict["Tin"]
    Tout = Tdict["Tout"]

    numlayers = len(Mlist)

    for m in range(numlayers):
        Mlistnewtemp = Mlist[numlayers - 1 - m]
        Tintemp = float(1.00)
        if (numlayers - 1 - m) == 0:
            pass
        else:
            for n in range(numlayers - 1, numlayers - m - 1, -1):
                if m == 0:
                    pass
                else:
                    for i in range(len(Tin[m])):
                        Ttemp = Tin[n][i]
                        Tintemp = Tintemp * Ttemp * Ttemp  # Note the square
            Mlistnewtemp = Mlistnewtemp * Tintemp
        Mlistnew1.append(Mlistnewtemp)

    Mlistnew1.reverse()

    for m in range(numlayers):
        Mlistnewtemp = Mlistnew1[m]
        Touttemp = float(1.00)
        for n in range(m, numlayers):
            Touttemp = Tout[n] * Touttemp
        Mlistnewtemp = Mlistnewtemp * Touttemp  # Note that this is NOT squared.
        Mlistnew.append(Mlistnewtemp)

    return Mlistnew
