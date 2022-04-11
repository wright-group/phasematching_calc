.. Lasers:

Lasers
=============

The :class:`Lasers` class consist of the needed input and output laser fields used in the four-wave mixing study,
as well as the methods for modifying these fields.  

Fields include the laser frequencies, their angles, their polarizations, the laser output combination, a standard
geometry the lasers make entering the sample, and some values helpful for computing phasematching effects 
generally related to the geometry.

The user builds on the Laser by modifying frequencies, angles, polarizations and output combination as lists,
then can change individual elements.   The user can change the input geometry, and can save or load the Lasers object.


Public methods are:

`addfrequencies(self, ar1)`
-------------------------
Input a list `ar1` of frequencies for each input.  As the supported calculations are meant for four-wave mixing,
the method only allows for three-member arrays/lists.  Input 1 is the first index of the list, 2 is the second,
3 is the third.

`changefrequency(self, freqnum, value)`
-------------------------
Changes the frequency of `freqnum` (input) to the `value`.  Should be a non-zero, real number.

`addangles(self, ar1)`
-------------------------
Input a list `ar1` of Cartesian angles in degrees in air entering the sample.  As the supported calculations are
meant for four-wave mixing, the method only allows for three-member arrays/lists.  Input 1 is the first index of the list, 2 is the second,
3 is the third.  The inputs are Cartesian...i.e., not Eulerian, because the expected geometries are to have
a laser input aligned along x or y relative to the sample focus but not both.   This simplifies Fresnel coefficient
calculations.

`changeangles(self, freqnum, value)`
-------------------------
Changes the anglee (in degrees) of `freqnum` (input) to the `value`.  Should be a real value below 90 degrees.

`addkcoeffs(self,ar1)`
----------------------
Input a list `ar1` of the four-wave mixing coefficient for each input.  Due to the limited number of 
interactions allowed in four-wave mixing, the expected set ranges integers [-3,3] including 0. However the
method does not restrict higher numbers due to certain allowances for the four-wave mixing program to 
effectively simulate certain higher-order ones.  No single element change method was created for this entry.

`addpolarizations(self,ar1)`
-----------------------------
Add an array of polarizations `ar1` with length equal to the input frequencies.
Value:  1 == Vertical,   !1 == Horizontal.   No other polarizations supported, and must be int.
No single element change method was created for this entry.   (Note: the output polarization is
"guessed" based on isotropic averaging)

`changegeometry(self,newval="boxcars")`
--------------------------------------
the `supportedgeometrylist` contains entries that this method may use.  The two important ones are a fully planar
or "2D" input/output geometry, or a "boxcars" (actually more like a "plus-sign") geometry, as they each
can be readily used for calculating Fresnel coefficients and each have their uses.  
        "boxcars":                           "planar":                                  Cartesian Axes:

                  O 2                           
                  |                                                                              |  y
         4        |         1                                                                    |
         O--------+---------O                ---------+----O-----O--O-------O                    +------ x
                  |                                        2     3  4       1
                  |
                  O 3  

        where "+" is the center of focus of the beams, and 4 is the output location. 

It is expected for reasonable thickness and angle calculations that the "+" represents the surface normal
to the IsoSample object used in tandem with this object.   

`save(self, file)`
------------------
Save the Lasers object to a JSON representation on disk.

`load(self, file)`
-----------------
Loads in the JSON disk file into the current object, modifying its values to that in the JSON file.


Other methods use masks to tell the phasematching calculator to skip computation of other Cartesian
angles as the inputs are aligned along specific coordinates.  As these above methods are used, another
methods reconverts the angles into radians for use in the phasematching calculation.




