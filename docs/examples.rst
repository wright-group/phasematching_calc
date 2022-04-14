.. examples:

Some examples using the ``phasematching_calc`` module are shown. 

Example 1. A 2D calculation of M factors and conversion into a ``WrightTools.data`` object.

First, the files are loaded into an IsoSample object:

.. plot::
    lay1file=os.path.join(filepath, 'CaF2_Malitson.txt')
    lay2file=os.path.join(filepath, 'H2O_1.txt')

    tkcaf2=0.02 
    tkwater=0.01 

    samp1=pc.IsoSample.IsoSample()
    desc="FWM cell"
    samp1.description=desc
    samp1.load_layer(lay1file, tkcaf2, label="caf2fw")
    samp1.load_layer(lay2file, tkwater, label="water")
    samp1.load_layer(lay1file, tkcaf2, label="caf2bw")

Then, the Lasers object is created from the Python commands (as opposed to being loaded from file):

.. plot::
    las=pc.Lasers.Lasers()
    arr1=[1800.0,2700.0,30000.0]
    las.add_frequencies(arr1)
    arr2=[20.0,7.0, 0.0]
    las.add_angles(arr2)
    arr3=[-1,1,1]
    las.add_k_coeffs(arr3)
    arr4=[1,1,1]
    las.add_polarizations(arr4)
    las.change_geometry("planar")

The above two objects signify we are looking at a DOVE process where the Lasers are at 1800, 2700, and 30000 cm-1,
the geometry is planar, and the first laser is the -k2 one while the others are k1 and k3.  The sample is a 
sandwich of caf2, water, and caf2 at thicknesses shown.

The simulated ``WrightTools`` Data conversion is shown in the next block.  One defines the numpy linspaces for the two independent
coordinates, the k1 and -k2 inputs.  Then a for loop cycles through each series of m,n elements and places it 
into a channel array that gets placed into a Data object.  The remaining code follows the WrightTools methodology
of using the linspaces as variables and plotting the result using the ``WrightTools.artists.quick2D`` method.


.. plot::
    var1=np.linspace(2450.00,2900.00,46)[:,None]
    var2=np.linspace(1300.0,1900.0,61)[None, :]
    var2a=np.linspace(1300.0,1900.0,61)

    ch1= np.zeros([len(var1), len(var2a)])
    for m in range(len(var1)):
        for n in range(len(var2a)):
            las.changefreq(1,var1[m])
            las.changefreq(2,var2a[n])
            Mlist,tklist,Tlist=pc.phasematch.m_calc(samp1,las)
            ch1[m,n]=np.abs(Mlist[1])


    data=wt.Data(name="FWM cell water w/CaF2 planar DOVE")
    data.create_variable(name="w1", units="wn", values= var1)
    data.create_variable(name="w2", units="wn", values= var2)
    data.create_channel(name='Mfactor', values=ch1)
    data.transform("w1","w2")
    wt.artists.quick2D(data)
    plt.show()


.. image:: Figure_1.png


**Example 2**. A similar calculation with a single 300 micron CaF2 window and different input geometries,
in this case a `boxcars` geometry.

.. plot::
    lay1file=os.path.join(filepath, 'CaF2_Malitson.txt')
    
    tkcaf2=0.03 #cm

    samp1=pc.IsoSample.IsoSample()
    desc="caf2window300um"
    samp1.description=desc
    samp1.load_layer(lay1file, tkcaf2, label="caf2")

    las=pc.Lasers.Lasers()
    arr1=[1800.0,2700.0,18400.0]
    las.add_frequencies(arr1)
    arr2=[8.0,8.0, 8.0]
    las.add_angles(arr2)
    arr3=[-1,1,1]
    las.add_k_coeffs(arr3)
    arr4=[1,1,1]
    las.add_pols(arr4)
    las.change_geometry("boxcars")

    var1=np.linspace(2600.00,3200.00,61)[:,None]
    var2=np.linspace(1600.0,2200.0,61)[None, :]
    var2a=np.linspace(1600.0,2200.0,61)

    ch1= np.zeros([len(var1), len(var2a)])
    for m in range(len(var1)):
        for n in range(len(var2a)):
            las.changefreq(1,var1[m])
            las.changefreq(2,var2a[n])
            Mlist,tklist,Tlist=pc.phasematch.m_calc(samp1,las)
            ch1[m,n]=np.abs(Mlist[0])  

    data=wt.Data(name="CaF2 300 micron boxcars DOVE")
    data.create_variable(name="w1", units="wn", values= var1)
    data.create_variable(name="w2", units="wn", values= var2)
    data.create_channel(name='Mfactor', values=ch1)
    data.transform("w2","w1")
    wt.artists.quick2D(data)
    plt.show()


.. image:: Figure_2.png

**Example 3**.  An angle solving routine for an oriented sapphire:acetonitrile:sapphire sample.
We assume the oriented sapphire limits its anisotropy to very small amounts that are neglected
and may approximate an isotropic sample.  This is reverting back to a planar geometry.  The Sympy
syntax requires the conversion of the `FiniteSet` to a `list`.  

The script is somewhat lengthy because finding angles is time-consuming and it is better to
use the "isclose" feature when possible.  One would need to calculate the angles crudely with
``isclose`` set to ``False`` then use the values calculated in further sections.  As the
``isclose`` feature only provides one solution, both angles solves would have to be separated.
The double for loops shown can be consolidated if one sets up separate ``Lasers`` objects.


.. plot::
    lay1file=os.path.join(filepath, 'CH3CN_paste_1.txt')
    lay2file=os.path.join(filepath, 'sapphire1.txt')
    lay3file=os.path.join(filepath, 'CaF2_Malitson.txt')

    tksap=0.02
    tkacn=0.01 
    tkcaf2=0.02

    # generation of a IsoSample
    samp1=pc.IsoSample.IsoSample()
    desc="FWM cell"
    samp1.description=desc
    #samp1.load_layer(lay1file, tksap, label="sapphire")
    samp1.load_layer(lay3file, tkcaf2, label="caf2")
    samp1.load_layer(lay1file, tkacn, label="acn")
    samp1.load_layer(lay3file, tkcaf2, label="caf2bw")
    #samp1.load_layer(lay1file, tksap, label="sapphire")

    #generation of a Lasers object.
    las=pc.Lasers.Lasers()
    arr1=[2200.0, 3150.0,17200.0]
    las.add_frequencies(arr1)
    arr2=[-13.0,6.0, 0.0]
    las.add_angles(arr2)
    arr3=[-1,1,1]
    las.add_k_coeffs(arr3)
    arr4=[1,1,1]
    las.add_pols(arr4)
    las.change_geometry("planar")

    var1=np.linspace(1600.0,2200.0,61)[None, :]
    var1a=np.linspace(1600.0,2200.0,61)
    var2=np.linspace(2600.00,3200.00,61)[:,None]
    var2a=np.linspace(1600.0,2200.0,61)

    ch1= np.zeros([len(var1a), len(var2a)])
    mold=int(0)
    for m in range(len(var1a)):
        for n in range(len(var2a)):
            las.change_freq(1,var1a[n])
            las.change_freq(2,var2a[m])
            if ((m==0) & (n==0)):
                angleair2=list(pc.phasematch.solve_angle(samp1,las,2,1,isclose=False))
                angletemp=angleair2[0]   # this needs to solve for remainder to work
                if np.any(angleair2):
                    ch1[m,n]=(angleair2)[0]
                    las.change_angle(1,angleair2[0])  
            elif (mold==m):
                angleair2=list(pc.phasematch.solve_angle(samp1,las,2,1,isclose=True))
                if np.any(angleair2):
                    ch1[m,n]=(angleair2)[0] 
                    las.change_angle(1,angleair2[0])           
            else:
                las.change_angle(1,angletemp) 
                angleair2=list(pc.phasematch.solve_angle(samp1,las,2,1,isclose=True))
                mold=m
                if np.any(angleair2):
                    ch1[m,n]=angleair2[0]
                    angletemp=angleair2[0]
                    las.change_angle(1,angleair2[0])

    data=wt.Data(name="angle for lower frequency input, opp side")
    data.create_variable(name="w1", units="wn", values= var1)
    data.create_variable(name="w2", units="wn", values= var2)
    data.create_channel(name='angleforw1', values=ch1)
    data.transform("w2","w1")
    wt.artists.quick2D(data)
    plt.show()


    for m in range(len(var1a)):
        for n in range(len(var2a)):
            las.change_freq(1,var1a[n])
            las.change_freq(2,var2a[m])
            if ((m==0) & (n==0)):
                angleair2=list(pc.phasematch.solve_angle(samp1,las,2,1,isclose=False))
                angletemp=angleair2[1]   # this needs to solve for remainder to work
                if np.any(angleair2):
                    ch1[m,n]=(angleair2)[1]
                    las.change_angle(1,angleair2[1])  
            elif (mold==m):
                angleair2=list(pc.phasematch.solve_angle(samp1,las,2,1,isclose=True))
                if np.any(angleair2):
                    ch1[m,n]=(angleair2)[0] 
                    las.change_angle(1,angleair2[0])           
            else:
                las.change_angle(1,angletemp) 
                angleair2=list(pc.phasematch.solve_angle(samp1,las,2,1,isclose=True))
                mold=m
                if np.any(angleair2):
                    ch1[m,n]=angleair2[0]
                    angletemp=angleair2[0]
                    las.change_angle(1,angleair2[0]) 

    data2=wt.Data(name="angle for lower frequency beam, same side")
    data2.create_variable(name="w1", units="wn", values= var1)
    data2.create_variable(name="w2", units="wn", values= var2)
    data2.create_channel(name='angleforw1', values=ch1)
    data2.transform("w2","w1")
    wt.artists.quick2D(data2)
    plt.show()

.. image:: Figure_3.png

.. image:: Figure_3b.png

Note the check is for the -k2 beam (i.e., "w1") and it is looking for phasematching in the acetonitrile layer (layernum=2).

The solution for the phasematching on the same side puts the two beams (w1 and w2)  at nearly identical angles.   This may
be good for certain optics, but bad if one wants separate optics for each beam.


**Example 4**.  A frequency solving routine for an oriented sapphire:acetonitrile:sapphire sample.
The conditions are virtually identical to Example 3 except that a frequency solve for the high frequency
k3 beam is requested.  The code is not posted as it is nearly identical except for  replacing the
line ``angleair2=pc.phasematch.solve_angle(samp1,las,2,1)`` with ``angleair2=pc.phasematch.solve_frequency(samp1,las,2,3)``.

.. image:: Figure_4.png

The expected w3 colors range from 18000 cm-1 at right to almost 19000 cm-1 at left, suggesting a very large change of colors
required that may obviate the method or require some additional laser modification for assistance.


**Example 5**.  A delta t check of the inputs in a thick sample between two caf2 windows.  A thick (1 mm) sample of
acetonitrile is simulated instead.  This thickness tends to be the upper limit for our liquid phase samples, as
 geometrical interactions tend to limit thicknesses.  (Geometrical calculations may be instituted as a function in a later version.)  

The code starts normally:

.. plot::
    lay3file=os.path.join(filepath, 'CaF2_Malitson.txt')
    lay4file=os.path.join(filepath, 'CH3CN_paste_1.txt')
    lay5file=os.path.join(filepath, 'CaF2_Malitson.txt')

    tkcaf2=0.02 
    tkacn=0.1 

    samp1=pc.IsoSample.IsoSample()
    desc="FWM cell"
    samp1.description=desc
    samp1.load_layer(lay5file, tkcaf2, label="caf2fw")
    samp1.load_layer(lay4file, tkacn, label="ACN")
    samp1.load_layer(lay3file, tkcaf2, label="caf2bw")

    las4=pc.Lasers.Lasers()
    arr1=[3150.0,2250.0,20000.0]
    las4.add_frequencies(arr1)
    arr2=[5.0,10.0,0.0]
    las4.add_angles(arr2)
    arr3=[1,-1,1]
    las4.add_k_coeffs(arr3)
    arr4=[1,1,1]
    las4.add_pols(arr4)
    las4.change_geometry("planar")

    tin,tout=pc.phasematch.calculate_ts(samp1,las4)
    print(tin,tout)

Some additonal code is needed to convert the times into more meaningful ones.  For example, the mean of
all 4 inputs and output was determined per layer, and the difference from that mean plotted per input.
.. plot::
    for m in range(len(tin)):
        if m == 0:
            pass
        else:
            for i in range(len(tin[m])):
                tin[m][i]=tin[m][i]-tin[m-1][i]

    for i in range(len(tout)):
        if i ==0:
            pass
        else:
            tout[i]=tout[i]-tout[i-1]

    print(tin,tout)
    tlist=list()
    x1=list()
    x2=list()
    x3=list()
    x4=list()
    y1=list()
    y2=list()
    y3=list()
    y4=list()

    for m in range(len(tin)):
        tinvec=list(tin[m])
        tinvec.append(tout[m])
        avg=np.mean(tinvec)
        tinvec=np.asarray(tinvec-avg)
        for i in range(len(tinvec)):
            if (i==0):
                x1.append(m+1)
                y1.append(tinvec[i])
            elif (i==1):
                x2.append(m+1)
                y2.append(tinvec[i])
            elif (i==2):
                x3.append(m+1)
                y3.append(tinvec[i])
            elif (i==3):
                x4.append(m+1)
                y4.append(tinvec[i])
            else:
                pass

    plt.rcParams['figure.autolayout']=True
    plt.xlim(0,5)
    plt.ylim(-60.0,30.0)
    plt.grid()

    xn1=x1
    yn1=y1
    plt.scatter(xn1,yn1, c="red")

    xn1=x2
    yn1=y2
    plt.scatter(xn1,yn1, c="green")

    xn1=x3
    yn1=y3
    plt.scatter(xn1,yn1, c="blue")

    xn1=x4
    yn1=y4
    plt.scatter(xn1,yn1, c="black")
    plt.show()


.. image:: Figure_5.png

Note the `scatter` plot does not show axes.  X is the layer number and y is the delta in femtoseconds each
input or output makes relative to the mean of the 4 at the end of the layer.  Red is input 1, green is input 2, 
blue is input 3, and black is the output.  At layer 3 input 1 and 2 or overlapped.   The first input makes a
40 fsec departure from the other inputs at the end of the acetonitrile layer, because the pulse is 
near a strong absorption that will "delay" it. However, in general all four are within 20 fsec of each other and
so with pulses 10x or wider in time there should be negligible effects on signal contributions due to delaying.

While a relative delay of 40 fsec is small for a frequency scanning method of four-wave mixing using picosecond or
similar pulses, it is occasionally useful to examine these delays, as often delays are established between pulses
to limit background signal.   The delays may be a simple pulsewidth difference.  In this case, ther could be small
changes in the expected delay during a scan by nature of the change in refractive index, which may slighly reduce
or increase the delay there.  At the edge of a pulsewidth, profound changes in background can occur, so noticeable
changes in background contributions may manifest in these areas.


**Example 6**.  A simple angle and frequency check.   Reverting back to the thin caf2:acetonitrile:caf2 sample,
a set of two frequency and angle solves are made for what may be considered two nearby data points to see 
how much of either should be made to achieve phasematching for both points.


.. plot::

filepath=os.path.join(ROOT_DIR, 'tests')

lay3file=os.path.join(filepath, 'CaF2_Malitson.txt')
lay4file=os.path.join(filepath, 'CH3CN_paste_1.txt')

tksapph=0.02 #cm
tkacn=0.01 #cm

samp1=pc.IsoSample.IsoSample()
desc="FWM cell"
samp1.description=desc
samp1.load_layer(lay3file, tksapph, label="caf2fw")
samp1.load_layer(lay4file, tkacn, label="ACN")
samp1.load_layer(lay3file, tksapph, label="caf2bw")

las4=pc.Lasers.Lasers()
arr1=[3150.0,2200.0,17200.0]
las4.add_frequencies(arr1)
arr2=[6.0,-13.20,0.0]
las4.add_angles(arr2)
arr3=[1,-1,1]
las4.add_k_coeffs(arr3)
arr4=[1,1,1]
las4.add_pols(arr4)
las4.change_geometry("planar")

angl1=pc.phasematch.solve_angle(samp1,las4,2,2, isclose=False)
out=list(angl1)
print(out[0])

freq=pc.phasematch.solve_frequency(samp1,las4,2,3,20)
out=list(freq)
print(out[0])

las4.change_freq(3,out[0])

las4.change_freq(2,2190.0)
angle=pc.phasematch.solve_frequency(samp1,las4,2,3,20)
out2=list(angle)
print(out2[0])

las4.change_freq(3,out[0])
angle=pc.phasematch.solve_angle(samp1,las4,2,2, isclose=False)
out3=list(angle)
print(out3[0])


Results are:
.. code-block:: python
 -13.2000000000000
17200.0000000000
17360.0000000000
-13.0000000000000

    
In this example, changing w3 by +160 cm-1 would result in the same phasematching as an angle change of +0.20 degrees.
Changes in w3 in this range would result in very large wavelength changes needed over an entire scan.  On the 
other hand, phasematching angle changes may be restricted to a small range due to aberrations.  It is possible
that the two can be modified in tandem in some studies.


**Example 7**.  Comparison of DOVE vs TSF signal intensity.  WIth the oriented sapphire:water:sapphire sample,
a check was done between the two expected signal intensities generated by the water layer in two example
four-wave mixing modes (DOVE vs TSF).  The H2O signal was not phasematcheable in DOVE with the w3 wavelength.  
However, it is important to note that as w3 increases, the vector contributions of k1 and -k2 become
smaller relative to k3, and so phasemismatching becomes less problematic for DOVE.  


.. plot::
    lay1file=os.path.join(filepath, 'sapphire1.txt')
    lay2file=os.path.join(filepath, "H2O_1.txt")
    tksap=0.02 
    tkwat=0.01

    samp1=pc.IsoSample.IsoSample()
    desc="sapphwatersapph"
    samp1.description=desc
    samp1.load_layer(lay1file, tksap, label="saphfw")
    samp1.load_layer(lay2file, tkwat, label="h2o")
    samp1.load_layer(lay1file, tksap, label="saphfw")

    las=pc.Lasers.Lasers()
    arr1=[1800.0,2700.0,30000.0]
    las.add_frequencies(arr1)
    arr2=[-18.0,8.0, 0.0]
    las.add_angles(arr2)
    arr3=[-1,1,1]
    las.add_k_coeffs(arr3)
    arr4=[1,1,1]
    las.add_pols(arr4)
    las.change_geometry("planar")

    var1=np.linspace(2450.00,2900.00,91)[:,None]
    var2=np.linspace(1300.0,1900.0,161)[None, :]
    var2a=np.linspace(1300.0,1900.0,161)

    ch1= np.zeros([len(var1), len(var2a)])
    ch2=np.zeros([len(var1), len(var2a)])
    ch3=np.zeros([len(var1), len(var2a)])

    for m in range(len(var1)):
        for n in range(len(var2a)):
            las.changefreq(1,var1[m])
            las.changefreq(2,var2a[n])
            Mlist,tklist,Tdict=pc.phasematch.m_calc(samp1,las)
            Alist, Alistout=pc.phasematch.calculate_absorbances(samp1,las) 
            Mlist1a=pc.phasematch.apply_absorbances(Mlist,Alist,Alistout)
            Mlist1b=pc.phasematch.apply_trans(Mlist1a, Tdict)
            ch1[m,n]=Mlist[1]  

    vec2=[1,1,1]
    las.addkcoeffs(vec2)

    for m in range(len(var1)):
        for n in range(len(var2a)):
            las.changefreq(1,var1[m])
            las.changefreq(2,var2a[n])
            Mlist2,tklist2,Tlist2=pc.phasematch.m_calc(samp1,las)
            ch2[m,n]=Mlist2[1]  

    ch3=ch1/ch2

    data=wt.Data(name="example")
    data.create_variable(name="w1", units="wn", values= var1)
    data.create_variable(name="w2", units="wn", values= var2)
    data.create_channel(name='DOVE', values=ch1)
    data.create_channel(name='TSF', values=ch2)
    data.create_channel(name="DOVE_TSF_RATIO", values=ch3)
    data.transform("w1","w2")
    wt.artists.quick2D(data, channel=0)
    plt.show()

    wt.artists.quick2D(data, channel=1)
    plt.show()

    wt.artists.quick2D(data, channel=2)
    plt.show()


.. image:: Figure_7a.png

.. image:: Figure_7b.png

.. image:: Figure_7c.png


Note the M factor is still quite large.  This calculation shows the effects of absorption within the water layer. 
The application of absorbances from later layers was not shown (they were calculated but not put
into the graphic.)

Here, the DOVE to TSF ratio can be up to a factor of 100 for this sample.  This is indicative of the expected
signal differences between the two processes strictly due to phase mismatching and not infrared or Raman
polarizabilities of compounds within the scan range.  A thin film of material would likely want to be added as 
an extra Layer, and the ratios between the two at that thin layer should approach 1 as it becomes small.

This comparison can be made with other samples.  Thicker, more transparent samples can yield DOVE/TSF ratios
into the 10^6 range.   