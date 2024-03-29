.. phasematching:

Phasematching
=============

The phasematching module consists of the functions used to calculate four-wave mixing phasematching effects in
multiplayer, isotropic samples.   The functions primarily require use of the :class:`phasematching_calc.IsoSample`
and :class:`phasmatching_calc.Lasers` objects, and so the user must have generated these objects prior.  Some
functions do not require these objects.

The functions described are meant for a "low-order" simulation of effects due to phase(mis)matching, absorption,
and time differences between inputs.   Higher-order phenomena such as non-linear refractive index, group velocity
changes (pulse chirping), Fresnel effects when layers become much smaller than wavelengths of light, etc. are not
incorporated.

Public functions are:


`m_calc (IsoSample, Lasers)`
-------------------------

Given the two objects, this function outputs the associated M factors as per equation 18b, Carlson and Wright,
_Applied Spectroscopy_ 43(7) 1195.  The M factors are calculated in each layer and are dependent on previous
layers by the Snell angle of refractions.   The M factor is a "squared norm" term and is therefore a real number.
The complex term found inside the term prior to squaring implies that each layer could have a "phase factor"
(exp(i*Delta)) relative to the previous layer(s), should more accurate modeling of phase effects between layers
be necessary.

Mcalc returns a tuple (1,2,3,4) of 4 outputs:

(1) an m-layered list of M factors for the specific phasematching path defined in the Lasers object.
See the examples in the scripts folders for further information.
(2) A phase for the term inside the M^2 expression, which can be incorporated into the square
root of the m-layered list to form a complex term, this phase is relative to the layer preceeding
it.
(3) a series of thicknesses experienced by the output four-wave mixing term, per layer.  The thickness
is the original normal thickness of the sample layer divided by cosine(launchangle).   This effective
thickness can be used to calculate four-wave mixing intensities more accurately.
(4) A dictionary of transmission coefficients based on Fresnel coefficients between each layer.  The key 'Tin'
represents the transmission coefficients for the inputs.  The first layer represents transmission from the air
into the sample, second represents transmission from 1st layer into 2nd layer, and so on.  It does not include
transmission between the final layer and the output air.   The key 'Tout' represents the transmission of the
output into the next layer.  Its first element starts with the transmission from layer 1 into layer 2 and ends
with the final layer transmission back into air.   The final key 'launchangledeg' is the Cartesian angle of
the output beam back into the air.  It presumes the input angles are all measured in a single pure Cartesian
coordinate relative to some central point such as the surface normal of the sample.

The Fresnel coefficients can be applied to incorporate reflection losses per layer in the sample.  It becomes
important to track whether these coefficients should be squared (as if they were used in the M factor which is
a "squared" term) or not (if one is just using the program to model beam transmissions, for example).  The "Tout"
terms have a different tracking consideration.  More to follow in the "applyfresneltrans" function.


`solve_angle(Iso,Las,layernum,freqnum, frequency=None, amt=None, isclose=False)`
----------------------------------------------------
Given the two objects, the functions attempts to solve through brute force iteration the angle that the
input at freqnum needs to make in the layernum to achieve phasematching.   The iteration proceeds through
directional and smaller amounts until convergence is met by an internal tolerance.  If `frequncy` is specified,
the function will work with that frequency in place of the one found in the Lasers object.  The brute force
algorithm starts at the angle originally specified in the Lasers object.

`solve_angle` returns a tuple: Sympy Set of Angles, amount.  It could be a single item `FiniteSet`, a double item
`FiniteSet`, an `Interval`, or an `EmptySet` if no solution is found.  An `Interval`` indicates a full range of
angles are possible, which can happen if the process can be phasematcheable in that manner.  The endpoint of the
interval is determined by critical angles that may be found between layers.  The return is for the original angle
in air of that input in degrees.  Amount is a value that can be used as an argument in later calls.

The keyword arg `amt` if not None directs the solver to work in increments given by its value.  This can
speed up the solving process for ``isclose``.  Use the amount from the tuple returns if needed.

The brute force algorithm is greatly speeded if ``isclose`` is set to ``True``.  However, only one solution will
be found this way.


`solve_frequency(Iso, Las, layernum, freqnum, amt=None, isclose=False)`
----------------------------------------------------
Given the two objects, the functions attempts to solve through brute force iteration the frequency of the input
required to achieve phasematching for the freqnum in the layernum.   The extra parameter `amt` is the amount of change
the brute force algorithm may be requested to move per iteration.  If not specified, an internal algorithm
approximates this amount.   This number may be required for a user to achieve faster convergences due to the
often large number of iterations required to find convergence in this function.  As a rule, it ought to be
roughly 1% of the original frequency of the input.

Setting `isclose` to `True` is not equivalent to that in `solve_angle`.  In this application, the tolerance is
improved.

`solve_frequency` returns a tuple:  Sympy Set of Frequencies, amount.  It could be a single item `FiniteSet`, an `Interval`,
or an `EmptySet` if no solution is found.  An `Interval`` indicates a full range of angles are possible, which can happen if
the process can be phasematcheable in that manner.  The endpoint is usually set to Infinity (`oo`), though a future
modification may end up limiting the frequency to a critical angle that may be found between layers.  See `solve_angle` for
information on what `amount` means.


`calculate_absorbances(Iso, Las)`
-------------------------------
Given the two objects, the function calculates log10 absorbances each input and output may make in each layer
of the sample.  This absorbance incorporates the angles the lasers make and assumes all original angles
were specified relative to the sample's surface normal (and that all layers are perfectly parallel).

It returns a tuple of lists `(Ain, Aout)`:  `Ain` specify the log10 absorbances of each laser through each layer,
while `Aout` is a list of the absorbances at the output for each layer.   These all incorporate the angles they
make in the sample layers according to Snell's Law.

This function can be used to calculate absorbance losses from earlier layers and how this affects four-wave
mixing intensity in the succeeding layer.   It may also serve as an auxiliary function for absorbance modeling
without need for use in four-wave mixing expressions.  See `apply_absorbances`.


`calculate_ts(Iso, Las)`
-------------------------------
Given the two objects, the function calculates the times each input and output make in each layer
of the sample.  This incorporates the angles the lasers make and their respective refractive indexes.

It returns a tuple of lists `(tin, tout)`:  `tin` specify the times in femtoseconds of each input makes by the end
of each layer, while `tout` is a list of times traversed by the output.

This calculation may be used to verify the pulses overlap properly in a sample or layer of a sample.
Ordinarily the differences are small and negligible relative to the pulse widths.  However, as the
phasematching calculator generalizes angles to very large numbers, and generalizes to very thick samples,
it is possible to find instances where the differences may manifest into a sizeable number that could be
important.


`apply_absorbances(Mlist, Alist_in, Alist_out=None)`
--------------------------------------------------
An auxiliary function not requiring the two IsoSample and Lasers objects.   It uses the Mfactor list from a
previous function and the absorbance lists from `calculate_absorbances`. The function calculates intensity
losses from prior absorbances into the succeeding layer's M factor.  It therefore squares each absorbance
loss as the M factor is a squared term.

If `Alist_out` is specified, the M-squared output four-wave mixing in that layer is scaled by a SINGLE, NON-SQUARED
absorbance sum of the SUCCESSIVE layers.  This is because it is presumed that four-wave mixing signal in a layer is
not dependent on four-wave mixing signal generated in previous layers...i.e., that the signal is "weak" relative
to the inputs.   (In other forms of non-linear mixing the generated signal may be a sizable fraction of the inputs
and then "rob" or "contribute" to the inputs prior...this would require knowledge of the phase factors described
above in `m_calc`.)

The output is a modified list of M factors `Mout` taking into account the absorbances.


`apply_trans(Mlist, Tdict=None)`
-------------------------------------
An auxiliary function not requiring the two IsoSample and Lasers objects.   It uses the Mfactor list from a
previous function and the fresnel `Tdict` from `m_calc`. The function calculates intensity
losses from prior reflections into the succeeding layer's M factor.  It therefore squares each transmission
coefficient as the M factor is a squared term.   The M-squared output four-wave mixing in that layer is also
scaled by a SINGLE, NON-SQUARED Fresnel reflection loss from successive layers.  This is because it is presumed
that four-wave mixing signal in a layer is not dependent on four-wave mixing signal generated in previous layers
...i.e., that the signal is "weak" relative to the inputs.   (In other forms of non-linear mixing the generated signal
may be a sizable fraction of the inputs and then "rob" or "contribute" to the inputs prior...this would require
knowledge of the phase factors described above in `m_calc`.)

The output is a modified list of M factors `Mout` taking into account the Fresnel losses.

As a zero-order calculation, this equation does not consider internal, interative reflections like in a cavity
at this time.


