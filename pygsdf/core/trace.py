#!/usr/bin/env python

# -*- coding: utf-8 -*-
r"""
Seismic SAC file Station class
==============================

read, process, and output single seismic traces

Class Sac: single station SAC format trace

Coded by: omid.bagherpur@gmail.com; github: omid-b
Update: 2023-10-13 (WORK in Progress!)

"""

import os

import obspy
import numpy as np


from obspy.core import UTCDateTime
from obspy.io.sac.sactrace import SACTrace
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics.base import locations2degrees

All_Possible_SAC_HEADERS = [\
    "A", "ADJTM", "AZ", "B", "BAZ", "CMPAZ", "CMPINC", "DELTA",
    "DEPMAX", "DEPMEN", "DEPMIN", "DIST", "E", "EVDP", "EVEL",
    "EVLA", "EVLO", "F", "GCARC", "IBODY", "IDEP", "IEVREG",
    "IEVTYP", "IFTYPE", "IINST", "IMAGSRC", "IMAGTYP", "INTERNAL",
    "IQUAL", "ISTREG", "ISYNTH", "IZTYPE", "KA", "KCMPNM", "KDATRD",
    "KEVNM", "KF", "KHOLE", "KINST", "KNETWK", "KO", "KSTNM", "KT0",
    "KT1", "KT2", "KT3", "KT4", "KT5", "KT6", "KT7", "KT8", "KT9",
    "KUSER0", "KUSER1", "KUSER2", "LCALDA", "LEVEN", "LOVROK", "LPSPOL",
    "MAG", "NEVID", "NORID", "NPTS", "NSNPTS", "NVHDR", "NWFID", "NXSIZE",
    "NYSIZE", "NZHOUR", "NZJDAY", "NZMIN", "NZMSEC", "NZSEC", "NZYEAR",
    "O", "ODELTA", "RESP0", "RESP1", "RESP2", "RESP3", "RESP4", "RESP5",
    "RESP6", "RESP7", "RESP8", "RESP9", "SB", "SDELTA", "STDP",
    "STEL", "STLA", "STLO", "T0", "T1", "T2", "T3", "T4", "T5",
    "T6", "T7", "T8", "T9", "USER0", "USER1", "USER2", "USER3",
    "USER4", "USER5", "USER6", "USER7", "USER8", "USER9", "XMAXIMUM",
    "XMINIMUM", "YMAXIMUM", "YMINIMUM"\
]


class Sac:
    """
    Class to read and process seismic data (SAC format; single SAC file)

    There are three main elements that could be accessed
    through their pertinent getters
        1) data >> get_data()
        2) times >> get_times()
        3) headers >> get_headers()

    """

    def __init__(self, sacfile=None,\
                 is_event=True, is_xcorr=False, headonly=False,):
        self.data = np.ndarray([])
        self.times = np.ndarray([])
        self.headers = {}
        self.starttime = None
        
        self.sacfile = sacfile
        self.headonly = headonly
        self.is_xcorr = is_xcorr
        if is_xcorr:
            self.is_event = False
        else:
            self.is_event = is_event
        
        if self.sacfile != None:
            self.read(self.sacfile)



    def read(self, sacfile):
        """
        Read sac file: append data, times, and headers
        """
        if os.path.isfile(sacfile):
            self.sacfile = sacfile
        else:
            self.sacfile = None
            print(f"Error: sac file does not exist: '{sacfile}'")
            exit(-1)

        # try read sac data and times
        try:
            st = obspy.read(self.sacfile, format="SAC", headonly=self.headonly)
            tr = st[0]
            if self.headonly:
                data = np.ndarray([])
                times = np.ndarray([])
            else:
                data = tr.data
                times = tr.times()
        except Exception as e:
            print(e)
            print(f"Error: could not read sac data: '{self.sacfile}'")
            exit(-1)

        headers = {}
        for key in  tr.stats.sac.keys():
            headers[f"{key}"] = tr.stats.sac[f"{key}"]

        self.append(data, times, headers)


    def append(self, data, times, headers):

        if type(headers) != dict:
            print("Error: appended 'headers' must be a type 'dict'.")
            exit(1)

        if len(data) ==  len(times) == 0:
            self.headonly = True
        else:
            self.headonly = False
            if len(times) != len(data):
                print("Error: number of elements in 'times' and 'data' do not match!")
                exit(1)

        self.data = data
        self.times = times
        self.headers = headers

        self.update_headers()



    def update_headers(self):
        
        required_headers = [
            'knetwk', 'kstnm', 'kcmpnm', 'stla', 'stlo', 'stel'
        ]
        for key in required_headers:
            if key not in self.headers.keys():
                print(f"Error: required header information is missing: '{key}'")
                exit(1)
        self.headers['knetwk'] = str(self.headers['knetwk'])
        self.headers['kstnm'] = str(self.headers['kstnm'])
        self.headers['kcmpnm'] = str(self.headers['kcmpnm'])
        self.headers['stla'] = float(self.headers['stla'])
        self.headers['stlo'] = float(self.headers['stlo'])
        self.headers['stel'] = float(self.headers['stel'])

        if not self.headonly:
            self.headers['delta'] = self.times[1] - self.times[0]
            self.headers['depmin'] = np.nanmin(self.data)
            self.headers['depmax'] = np.nanmax(self.data)
            self.headers['b'] = np.min(self.times)
            self.headers['e'] = np.max(self.times)
            self.headers['npts'] = len(self.times)

        starttime_headers = [
            'nzyear', 'nzjday', 'nzhour', 'nzmin', 'nzsec', 'nzmsec'
        ]
        if 'starttime' in self.headers.keys() and self.starttime == None:
            if type(self.headers['starttime']) == str:
                self.starttime = UTCDateTime(self.headers['starttime'])
            elif type(self.headers['starttime']) == UTCDateTime:
                self.starttime = self.headers['starttime']
            else:
                self.starttime = UTCDateTime("1970-01-01T00:00:00.000")
            del self.headers['starttime']

        elif all(x in list(self.headers.keys()) for x in starttime_headers):
            
            self.starttime = UTCDateTime(
                year = int(self.headers['nzyear']),
                julday = int(self.headers['nzjday']),
                hour = int(self.headers['nzhour']),
                minute = int(self.headers['nzmin']),
                second = int(f"{self.headers['nzsec']}"),
                microsecond = int(self.headers['nzmsec'] * 1000)
            )


        else:
            
            self.starttime = UTCDateTime("1970-01-01T00:00:00.000")

        self.headers['nzyear'] = int(self.starttime.year)
        self.headers['nzjday'] = int(self.starttime.julday)
        self.headers['nzhour'] = int(self.starttime.hour)
        self.headers['nzmin'] = int(self.starttime.minute)
        self.headers['nzsec'] = int(self.starttime.second)
        self.headers['nzmsec'] = int(self.starttime.microsecond / 1e3)

        # self.headers['iftype'] = int(1)
        del self.headers['iftype']
        self.headers['leven'] = int(1)
        self.headers['lpspol'] = int(0)
        self.headers['lovrok'] = int(1)
        self.headers['lcalda'] = int(1)

        # set event headers (?)
        if  'evla' in self.headers.keys()\
        and 'evlo' in self.headers.keys():
            dist, az, baz = gps2dist_azimuth(self.headers['evla'], self.headers['evlo'],\
                                             self.headers['stla'], self.headers['stlo'])
            gcarc = locations2degrees(self.headers['evla'], self.headers['evlo'],\
                                             self.headers['stla'], self.headers['stlo'])
            self.headers['az'] = float(az)
            self.headers['baz'] = float(baz)
            self.headers['dist'] = float(dist / 1000)
            self.headers['gcarc'] = float(gcarc)
            if 'evdp' not in self.headers.keys():
                self.headers['evdp'] = float(0.0)

            if 'o' not in self.headers.keys():
                self.headers['evdp'] = float(0.0)

        for key in self.headers.keys():
            if key.upper() not in All_Possible_SAC_HEADERS:
                del self.headers[key]


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
        str_headers = f"headers = %s{str_headers[0:-2]}%s" %('{', '}')

        return f"{str_data}\n{str_times}\n{str_headers}"


    def write(self, sacfile):
        if os.path.splitext(sacfile)[1] != '.sac':
            sacfile = f"{sacfile}.sac"
        tr = SACTrace( data=np.array(self.data), **self.headers)
        tr.write(sacfile)


    def __repr__(self):
        sac_object = {
            "data": self.data,
            "times": self.times,
            "headers": self.headers,
            "starttime": self.starttime
        }
        return sac_object


    def get_data(self):
        return self.data


    def get_times(self):
        return self.times


    def get_headers(self):
        return self.headers

    def get_starttime(self):
        return self.starttime

