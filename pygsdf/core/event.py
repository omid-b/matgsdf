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
            "otime":     <float: event_origin_time>,   # event origin time (float value in time axis)
            "otime_str": <str: event_origin_time>,     # event origin time in string in YYJJJHHMMSS format
            "otime_utc": <obspy UTCDateTime: event_origin_time>, # event origin time obspy UTCDateTime object
            "stations": <list(Station object)>, each station object has three main attributes: times, data, headers
        }

    """
    def __init__(self, event_dir=None, extensions=["sac", "SAC"]):
        self.evla = None
        self.evlo = None
        self.evdp = None
        self.otime = None
        self.otime_str = None
        self.otime_utc = None
        self.stations = []
        self.extensions = extensions
        self.read_folder(event_dir=event_dir,\
                         extensions=extensions)


    def read_folder(self, event_dir, extensions=None):
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
            self.append(event_sacfiles)



    def append(self, sacfiles):
        if type(sacfiles) == str:
            sacfiles = [sacfiles]
        
        for sacfile in sacfiles:
            sta = Station(sacfile)
            if not len(self.stations):
                self.evla = sta.get_headers()["evla"]
                self.evlo = sta.get_headers()["evlo"]
                self.evdp = sta.get_headers()["evdp"]
                self.utc_ = sta.get_headers()["o"]

                self.otime_str = sta.get_headers()["evla"]
                self.otime_utc = sta.get_headers()["evla"]
                self.stations = []





        



