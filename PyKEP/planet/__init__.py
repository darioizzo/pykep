# -*- coding: iso-8859-1 -*-

"""
This module contains all the classes that allow to construct efficiently
low-thrust tajectories using our own flavour of the Sims-Flanagan model: a trajectory
transcription method that forms the basis for MALTO, the software in use in JPL
for preliminary interplanetary trajectory design.
"""
from PyKEP.planet._planet import *
from PyKEP.planet._planet import _base
import PyKEP
import cPickle as pickle
if PyKEP.__extensions__["scipy"]:
    from scipy.interpolate import UnivariateSpline

def _keplerian_ctor(self, *args):
    """
PyKEP.planet.keplerian(when,orbital_elements, mu_central_body, mu_self,radius, safe_radius [, name = 'unknown'])

PyKEP.planet.keplerian(when,r,v, mu_central_body, mu_self,radius, safe_radius [, name = 'unknown'])

- when: a :py:class:`PyKEP.epoch` indicating the orbital elements epoch
- orbital_elements: a sequence of six containing a,e,i,W,w,M (SI units, i.e. meters and radiants)
- r,v: position and velocity of an object at when (SI units)
- mu_central_body: gravity parameter of the central body (SI units, i.e. m^2/s^3)
- mu_self: gravity parameter of the planet (SI units, i.e. m^2/s^3)
- radius: body radius (SI units, i.e. meters)
- safe_radius: mimimual radius that is safe during a fly-by of the planet (SI units, i.e. m)
- name: body name

.. note::

   All classes having Keplerian ephemerides as :py:class:`PyKEP.planet.mpcorb` inherit from this (c++) class

Example::

  earth = planet(epoch(54000,"mjd"),(9.99e-01 * AU, 1.67e-02, 8.85e-04 * DEG2RAD, 1.75e+02 * DEG2RAD, 2.87e+02 * DEG2RAD, 2.57e+02 * DEG2RAD), MU_SUN, 398600e9, 6378000, 6900000,  'Earth')"
    """
    self._orig_init(*args)
keplerian._orig_init = keplerian.__init__
keplerian.__init__ = _keplerian_ctor

def planet_spice(target=None, observer=None, ref_frame=None, aberrations=None, mu_central_body=None, mu_self=None, radius=None, self_radius=None):
    if not PyKEP.__extensions__["pygmo"]:
        print "PyGMO is required to use SPICE ephemerides, but not detected. Planet creation aborted."
        return None
    #if not spice_kernel_loaded():
        #print "There's no SPICE kernel loaded for instantiating a planet. Use PyKEP.util.load_spice_kernel() first."
        #return None

    missing = [arg for arg in locals().keys() if locals()[arg] is None]
    if missing:
        print "Missing parameters: %s." % ", ".join(missing)
        return None

    return PyKEP.planet.spice(target, observer, ref_frame, aberrations, mu_central_body, mu_self, radius, self_radius)keplerian

def dump(obj, fname="eph.p"):
    pickle.dump(obj, open(fname, "wb"))
def load(fname):
    return pickle.load(open(fname, "rb"))

class Interp():
    """
PyKEP.planet.Interp(planet, epochs [, samples=100, method = 'linear'])

An object that can be serialized with dump(obj) and load(fname).

-planet: The planet loaded from a SPICE kernel that will be serialized.
-epochs: A list of [start, end] to specify the range of epochs to save on serialization.
-samples: The amount of data points to extract from the epoch range.
-method: Interpolation method. Currently only 'linear' and 'spline' are supported, and 'spline' requires Scipy to be installed.

Example:: 

    ip = Interp(spice_pl, [17121, 19121.02], samples=100, method="spline")
    dump(ip, fname="eph.p")
    # Later:
    pl = load("eph.p")
    print pl.eph(17777.05)
    """
    def __init__(self, planet, epochs, samples=100, method="linear"):
        self.planet = planet
        self.epochs = epochs
        self.samples = samples  # Grid size for extracted epochs
        self.splines = []   # In the event of spline interpolation, we create splines using Scipy
        self.ser_eph = {}   # Serialized ephemeris data points between epochs
        self.indices = []   # The epochs corresponding to serialized ephemerides (for convenience)

        # Initialized on creation
        self.cache_eph()
        self.set_method(method)
        self.interpolate()

    def set_method(self, method):
        ''' Sets the interpolation method. '''
        if method=="spline":
            if not PyKEP.__extensions__["scipy"]:
                print "Scipy is required for spline interpolation. Reverting to default linear."
                method = "linear"
            
        elif method not in ["spline", "linear"]:
            print "Currently only spline and linear interpolations are supported. Reverting to default linear."
            method = "linear"

        self.method = method

    def cache_eph(self):
        ''' Samples ephemeris data points and the corresponding epochs in a grid for serialization. '''
        ds = (self.epochs[1]-self.epochs[0])/float(self.samples)

        # This range includes the start- and endpoints of epochs
        for s in range(self.samples+1):
            i = self.epochs[0]+ds*s
            self.ser_eph[i] = self.planet.eph(i)
            self.indices.append(i)

    def interpolate(self):
        ''' Readies the interpolation method for ephemeris requests. '''
        if self.method == "linear":
            # Linearly interpolated points are naively calculated upon request
            return
        elif self.method == "spline":
            # For this one, we create actual spline objects for later reference
            for sp in range(2):
                for sp2 in range(3):
                    y = [self.ser_eph[i][sp][sp2] for i in self.indices]
                    self.splines.append(UnivariateSpline(self.indices,y,s=1))
            
    def eph(self, eph):
        ''' Fetches the requested ephemeris through interpolation. '''
        if eph<self.epochs[0] or eph>self.epochs[1]:
            print "Requested ephemeris outside of the saved interpolation range."
            return None

        if eph in self.ser_eph:
            return self.ser_eph[eph]

        if self.method == "linear":
            for i,v in enumerate(self.indices):
                if v>eph:
                    next = v
                    prev = self.indices[i-1]
                    break

            # Eph location between data points normalized to 0..1, from which the interpolation
            dist = (eph-prev)/float(next-prev)

            # Faithful to the original format of the ephemerides 
            out = [[],[]]
            for d in range(2):
                for d2 in range(3):
                    prev_val = self.ser_eph[prev][d][d2]
                    next_val = self.ser_eph[next][d][d2]
                    interp = prev_val + dist * (next_val - prev_val)
                    out[d].append(interp)

            out = [tuple(o) for o in out]
            return out
            
        elif self.method == "spline":
            # Fetches these from saved splines
            # Faithful to the original format of the ephemerides 
            return [tuple([float(self.splines[s](eph)) for s in range(3)]),
                    tuple([float(self.splines[s2](eph)) for s2 in range(3,6)]) ]
 
