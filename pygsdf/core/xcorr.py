#!/usr/bin/env python

# -*- coding: utf-8 -*-
r"""
Seismic SAC file XCorr class
==============================
XCorr
Seismic SAC file XCorr class is to read, process cross-correlograms
e.g., Empirical Green's Function datasets. It is very similar to 
the Event class (in a sense that in contains station objects),
except that it also allows station objects with different
event locations. Moreover, the Event class allows for different
components in a single event object e.g. BHZ, HHZ, LHZ etc., hence 
one has to be very careful not appending unwanted sac files into a 
Event object. This limitation does not apply to the XCorr class as 
it only accepts exact same component read from header 'kcmpnm'.
Alternatively, if header 'kcmpnm' is not available one can specify 
parameter 'component' manually (at own risk!).

Coded by: omid.bagherpur@gmail.com; github: omid-b
Update: 2023-10-13 (WORK in Progress!)

"""

import os

import obspy
import numpy as np

from .sac import Sac

class XCorr:
    """
    Seismic SAC file XCorr class is to store Sac classes
    with the consistent XCorr information.

    XCorr (cross correlogram) object
    ------------
    An XCorr object is a python dictionary:
        XCorr = {
            "component": <float: XCorrs' component (either ZZ, RR, TT)>
            "tmin":      <float: minimum time>
            "tmax":      <float: maximum time>
            "delta":     <float: delta/sampling interval>
            "min_dist":  <float: minimum inter-station distance (km)> 
            "max_dist":  <float: maximum inter-station distance (km)> 
            "avg_dist":  <float: average inter-station distance (km)> 
            "sym":   <bool: all cross corellograms are symmetrized>
            "sacs": <list(Sac objects)>
        }

        component is read from header 'kcmpnm' and must be consistent throughout all sacs

    """
    def __init__(self, xcorr_dir=None, extensions=["sac", "SAC"],
                 sort_by_dist = True, component=None):
        self.component = component
        self.tmin = None
        self.tmax = None
        self.delta = None
        self.min_dist = None
        self.max_dist = None
        self.avg_dist = None
        self.sym = None
        self.sacs = []
        self.egf_object = {
            "component": self.component,
            "tmin": self.tmin,
            "tmax": self.tmax,
            "delta": self.delta,
            "min_dist": self.min_dist,
            "max_dist": self.max_dist,
            "avg_dist": self.avg_dist,
            "sym": self.sym,
            "sacs": self.sacs,
        }
        self.extensions = extensions
        self.read_folder(xcorr_dir=xcorr_dir,\
                         extensions=extensions,\
                         sort_by_dist=sort_by_dist)


    def read_folder(self, xcorr_dir, extensions=None, sort_by_dist=True):
        """
        Read SAC files within an cross correlogram dataset directory and append to self

        """
        if not extensions:
            extensions = self.extensions

        if not xcorr_dir:
            return
        elif not os.path.isdir(xcorr_dir):
            print(f"Error: input XCorr directory does not exist: '{xcorr_dir}'")
            exit(-1)

        xcorr_sacfiles = []
        for f in os.listdir(xcorr_dir):
            _, ext = os.path.splitext(f)
            ext = ext[1:]
            if ext in extensions:
                xcorr_sacfiles.append(os.path.join(xcorr_dir,f))

        if not len(xcorr_sacfiles):
            print(f"WARNING: No SAC file found in xcorr directory: '{xcorr_dir}'\nSAC file extensions: {' '.join(extensions)}")
            return
        else:
            self.append(xcorr_sacfiles, sort_by_dist=sort_by_dist)


    def append(self, sacfiles, sort_by_dist=True):
        if type(sacfiles) == str:
            sacfiles = [sacfiles]
        

        for sacfile in sacfiles:
            sac = Sac(sacfile, is_xcorr=True, component=self.component)
            sac.times += sac.headers['b'] # neccessary for cross correlation data!
            if not len(self.sacs):
                self.component = sac.headers['kcmpnm']
                self.tmin = sac.headers['b']
                self.tmax = sac.headers['e']
                self.delta = sac.headers['delta']
                self.sacs.append(sac)

            last_sac = self.sacs[-1]
            if  sac.headers['kcmpnm'] == last_sac.headers['kcmpnm'] and \
                sac.times[0] == last_sac.times[0] and \
                sac.times[-1] == last_sac.times[-1] and \
                sac.headers['delta'] == last_sac.headers['delta']:
                self.sacs.append(sac)
            
                
            else:
                print(f"Data append faild: '{sacfile}'")
                print(" >> Headers do not match for station '%s' and '%s'." \
                    %(sac.headers['kstnm'], last_sac.headers['kstnm']))
 
        if sort_by_dist:
            self.sort_by_distance()

        self.update_attributes()
        self.event_object = {
            "component": self.component,
            "tmin": self.tmin,
            "tmax": self.tmax,
            "delta": self.delta,
            "min_dist": self.min_dist,
            "max_dist": self.max_dist,
            "avg_dist": self.avg_dist,
            "sym": self.sym,
            "sacs": self.sacs,
        }


    def sort_by_distance(self):
        sacs_dists = []
        for ista in range(len(self.sacs)):
            sacs_dists.append(self.sacs[ista].headers["dist"])
        self.sacs = [x for _, x in sorted(zip(sacs_dists, self.sacs))]


    def update_attributes(self):
        # update distance and symmetry status attributes
        nsta = len(self.sacs)
        # update: min_dist, max_dist, avg_dist 
        dists = []
        for ista in range(nsta):
            dists.append(self.sacs[ista].headers['dist'])
        self.min_dist = np.nanmin(dists)
        self.max_dist = np.nanmax(dists)
        self.avg_dist = np.nanmean(dists)
        # update sym
        sym = True

        for ista in range(nsta):
            sac = self.sacs[ista]
            # decide based on the minimum time
            if self.tmin == 0 or self.tmax == 0:
                break # one-sided xcorr is always symmetric!
            # compare causal (positive lags) and acausal (negative lags) data
            causal = []
            acausal = []
            for itime, time in enumerate(sac.times):
                if time < 0:
                    acausal.append(sac.data[itime])
                elif time > 0:
                    causal.append(sac.data[itime])
            if len(causal) == 0 or len(acausal) == 0:
                break # one-sided xcorr is always symmetric!
            acausal = acausal[::-1] # reverse order
            if causal != acausal:
                sym = False
                break
        self.sym = sym


    def __str__(self):
        return "class object: <XCorr object>"


    def __repr__(self):
        return self.event_object