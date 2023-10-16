#!/usr/bin/env python

# -*- coding: utf-8 -*-
r"""
Seismic SAC file Event class
==============================

Seismic SAC file Event class is to read, process event data 
and contains Station objects with similar event header information

Coded by: omid.bagherpur@gmail.com; github: omid-b
Update: 2023-10-13 (WORK in Progress!)

"""

import os

import obspy
import numpy as np

from .station import Station


class Event:
    """
    Seismic SAC file Event class is to store Station classes
    with the same event information. It does not accept Station objects
    with inconsistent event information.

    Event object
    ------------
    An event object is a python dictionary:
        event = {
            "evla":  <float: event_latitude>,  # event latitude
            "evlo":  <float: event_longitude>, # event longitude
            "evdp":  <float: event_depth>, # event depth
            "otime":     <float: event_origin_time>,   # event origin time (SAC header 'o')
            "otime_str": <str: event_origin_time>,     # event origin time in string in YYJJJHHMMSS format
            "otime_utc": <obspy UTCDateTime: event_origin_time>, # event origin time obspy UTCDateTime object
            "stations": <list(Station object)>, each station object has three main attributes: times, data, headers
        }

    """
    def __init__(self, event_dir=None, extensions=["sac", "SAC"], sort_by_dist = True):
        self.evla = None
        self.evlo = None
        self.evdp = None
        self.delta = None
        self.starttime_utc = None
        self.endtime_utc = None
        self.otime_utc = None
        self.otime_str = None
        self.stations = []
        self.event_object = {
            "evla": self.evla,
            "evlo": self.evlo,
            "evdp": self.evdp,
            "delta": self.delta,
            "starttime_utc": self.starttime_utc,
            "endtime_utc": self.endtime_utc,
            "otime_utc": self.otime_utc,
            "otime_str": self.otime_str,
            "stations": self.stations,
        }
        self.extensions = extensions
        self.read_folder(event_dir=event_dir,\
                         extensions=extensions,\
                         sort_by_dist=sort_by_dist)


    def read_folder(self, event_dir, extensions=None, sort_by_dist=True):
        """
        Read/append SAC files within an event directory and add to self

        """
        if not extensions:
            extensions = self.extensions

        if not event_dir:
            return
        elif not os.path.isdir(event_dir):
            print(f"Error: input event directory does not exist: '{event_dir}'")
            exit(-1)

        event_sacfiles = []
        for f in os.listdir(event_dir):
            _, ext = os.path.splitext(f)
            ext = ext[1:]
            if ext in extensions:
                event_sacfiles.append(os.path.join(event_dir,f))

        if not len(event_sacfiles):
            print(f"WARNING: No SAC file found in event directory: '{event_dir}'\nSAC file extensions: {' '.join(extensions)}")
            return
        else:
            self.append(event_sacfiles, sort_by_dist=sort_by_dist)



    def append(self, sacfiles, sort_by_dist=True):
        if type(sacfiles) == str:
            sacfiles = [sacfiles]
        
        for sacfile in sacfiles:
            sta = Station(sacfile)
            if not len(self.stations):
                self.evla = sta.headers['evla']
                self.evlo = sta.headers['evlo']
                self.evdp = sta.headers['evdp']
                self.delta = sta.headers['delta']
                self.starttime_utc = sta.headers['starttime_utc']
                self.endtime_utc = sta.headers['endtime_utc']
                self.otime_utc = sta.headers['otime_utc']
                self.otime_str = sta.headers['otime_str']
                self.stations.append(sta)

            last_sta = self.stations[-1]
            if  sta.headers['evla'] == last_sta.headers['evla'] and \
                sta.headers['evlo'] == last_sta.headers['evlo'] and \
                sta.headers['evdp'] == last_sta.headers['evdp'] and \
                sta.headers['delta'] == last_sta.headers['delta'] and \
                sta.headers['starttime_utc'] == last_sta.headers['starttime_utc'] and \
                sta.headers['endtime_utc'] == last_sta.headers['endtime_utc'] and \
                sta.headers['otime_utc'] == last_sta.headers['otime_utc'] and \
                sta.headers['otime_str'] == last_sta.headers['otime_str']:
                self.stations.append(sta)
            else:
                print(f"Data append failed: '{sacfile}'")
                print(" >> Headers do not match for station '%s' and '%s'." \
                    %(sta.headers['kstnm'], last_sta.headers['kstnm']))
                continue

        if sort_by_dist:
            self.sort_by_distance()

        self.event_object = {
            "evla": self.evla,
            "evlo": self.evlo,
            "evdp": self.evdp,
            "delta": self.delta,
            "starttime_utc": self.starttime_utc,
            "endtime_utc": self.endtime_utc,
            "otime_utc": self.otime_utc,
            "otime_str": self.otime_str,
            "stations": self.stations,
        }

    def sort_by_distance(self):
        stations_dists = []
        for ista in range(len(self.stations)):
            stations_dists.append(self.stations[ista].headers["dist"])
        self.stations = [x for _, x in sorted(zip(stations_dists, self.stations))]

    
    def __str__(self):
        return self.otime_str


    def __repr__(self):
        return self.event_object

        






        



