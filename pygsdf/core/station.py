#!/usr/bin/env python

# -*- coding: utf-8 -*-
r"""
Seismic SAC file Station class
==============================

Seismic SAC file Station class is to read, process, and output 
single station SAC file data and headers

Coded by: omid.bagherpur@gmail.com; github: omid-b
Update: 2023-10-13 (WORK in Progress!)

"""

import os

import obspy
import numpy as np


class Station:
    """
    Class to read and process seismic data (SAC format)

    There are three main elements that could be accessed
    through their pertinent getters
        1) data >> get_data()
        2) times >> get_times()
        3) headers >> get_headers()

    """

    def __init__(self, sacfile=None, headers_only=False, \
                 event_headers=True, is_xcorr=False, component=None):
        self.data = None
        self.times = None
        self.headers = None
        self.filtered = None
        self.frequencies = None
        self.station_object = {
            "data": self.data,
            "times": self.times,
            "headers": self.headers,
            "filtered": self.filtered,
            "frequencies": self.frequencies,
        }
        self.sacfile = sacfile
        self.headers_only = headers_only
        self.is_xcorr = is_xcorr
        self.component = component
        if is_xcorr:
            self.event_headers = False
        else:
            self.event_headers = event_headers
        if self.sacfile:
            self.read(self.sacfile)


    def read(self, sacfile):
        """
        Read sac file data and headers
        """
        if os.path.isfile(sacfile):
            self.sacfile = sacfile
        else:
            self.sacfile = None
            print(f"Error: sac file does not exist: '{sacfile}'")
            exit(-1)

        # try read sac data and times
        try:
            st = obspy.read(self.sacfile, format="SAC", headonly=self.headers_only)
            tr = st[0]
            if self.headers_only:
                self.data = None
                self.times = None
            else:
                self.data = tr.data
                self.times = tr.times()
        except Exception as e:
            print(e)
            print(f"Error: could not read sac data: '{self.sacfile}'")
            exit(-1)

        # try reading sac headers
        try:
            knetwk = str(tr.stats.sac.knetwk)
            kstnm = str(tr.stats.sac.kstnm)
            if self.component:
                kcmpnm = self.component
            else:
                kcmpnm = str(tr.stats.sac.kcmpnm)
            stla = float(tr.stats.sac.stla)
            stlo = float(tr.stats.sac.stlo)
            try:
                stel = float(tr.stats.sac.stel)
            except:
                stel = -12345

            if self.event_headers:
                evla = float(tr.stats.sac.evla)
                evlo = float(tr.stats.sac.evlo)
                evdp = float(tr.stats.sac.evdp)
                gcarc = float(tr.stats.sac.gcarc)
                dist = float(tr.stats.sac.dist)
                az = float(tr.stats.sac.az)
                baz = float(tr.stats.sac.baz)
                o = float(tr.stats.sac.o)
                otime_utc = tr.stats.starttime + o
                otime_str = "%s%03d%02d%02d%02d" \
                            %(str(otime_utc.year)[2:],
                              otime_utc.julday,
                              otime_utc.hour,
                              otime_utc.minute,
                              otime_utc.second)
            elif self.is_xcorr:
                evla = float(tr.stats.sac.evla)
                evlo = float(tr.stats.sac.evlo)
                try:
                    evdp = float(tr.stats.sac.evdp)
                except:
                    evdp = -12345
                gcarc = float(tr.stats.sac.gcarc)
                dist = float(tr.stats.sac.dist)
                az = float(tr.stats.sac.az)
                baz = float(tr.stats.sac.baz)
                o = 0.
                otime_utc = tr.stats.starttime
                otime_str = "%s%03d%02d%02d%02d" \
                            %(str(otime_utc.year)[2:],
                              otime_utc.julday,
                              otime_utc.hour,
                              otime_utc.minute,
                              otime_utc.second)
            else:
                evla = None
                evlo = None
                evdp = None
                gcarc = None
                dist = None
                az = None
                baz = None
                o = None
                otime_utc = None
                otime_str = None

            b = float(tr.stats.sac.b)
            e = float(tr.stats.sac.e)
            delta = float(tr.stats.sac.delta)
            starttime_utc = tr.stats.starttime
            endtime_utc = tr.stats.endtime
            tag = f"{knetwk}.{kstnm}.{kcmpnm}"

            self.headers = {}
            self.headers['tag'] = tag
            self.headers['knetwk'] = knetwk
            self.headers['kstnm'] = kstnm
            self.headers['kcmpnm'] = kcmpnm
            self.headers['stla'] = stla
            self.headers['stlo'] = stlo
            self.headers['stel'] = stel
            self.headers['evla'] = evla
            self.headers['evlo'] = evlo
            self.headers['evdp'] = evdp
            self.headers['gcarc'] = gcarc
            self.headers['dist'] = dist
            self.headers['az'] = az
            self.headers['baz'] = baz
            self.headers['o'] = o
            self.headers['b'] = b
            self.headers['e'] = e
            self.headers['delta'] = delta
            self.headers['starttime_utc'] = starttime_utc
            self.headers['endtime_utc'] = endtime_utc
            self.headers['otime_utc'] = otime_utc
            self.headers['otime_str'] = otime_str
            self.headers['tag'] = tag
        except Exception as e:
            print(f"Error: could not read sac header '{e}': '{self.sacfile}'")
            exit(-1)

        self.station_object = {
            "data": self.data,
            "times": self.times,
            "headers": self.headers,
        }


    def __str__(self):
        str_data = ""
        for d in self.data:
            str_data += f"{str(d)}, "
        str_data = f"data = [{str_data[0:-2]}]"

        str_times = ""
        for t in self.times:
            str_times += f"{str(t)}, "
        str_times = f"times = [{str_times[0:-2]}]"

        str_headers = ""
        for hdr in self.headers.keys():
            if hdr in ['tag', 'knetwk', 'kstnm', 'kcmpnm']:
                str_headers += f'   "%s": "%s",\n' %(hdr, self.headers[f"{hdr}"])
            elif hdr == "utc":
                pass
            else:
                str_headers += f'   "%s": %s,\n' %(hdr, self.headers[f"{hdr}"])
        str_headers = f"headers = %s\n{str_headers[0:-2]}\n%s" %('{', '}')

        return f"{str_data}\n{str_times}\n{str_headers}"


    def __repr__(self):
        return self.station_object


    def get_data(self):
        return self.data


    def get_times(self):
        return self.times


    def get_headers(self):
        return self.headers

