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

import matplotlib.pyplot as plt
import numpy as np
import obspy
import scipy

from obspy.core import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics.base import locations2degrees
from obspy.io.sac.sactrace import SACTrace
from scipy.fft import fft, ifft, rfft, irfft

from scipy.signal import butter,filtfilt

All_Possible_SAC_Headers = [\
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

Required_SAC_Headers = [\
    'knetwk', 'kstnm', 'kcmpnm', 'stla', 'stlo', 'stel' \
]

class Trace:
    """
    Class to read and process seismic trace data 
    (continous timeseries; output file format: SAC)

    There are four main elements that could be accessed
    through their pertinent getter methods
        1) data >> get_data()
        2) times >> get_times()
        3) headers >> get_headers()
        3) starttime >> get_starttime()
        3) endtime >> get_endtime()

    """

    def __init__(self, sacfile=None):
        self.data = np.ndarray([])
        self.times = np.ndarray([])
        self.headers = {}
        self.starttime = None
        self.endtime = None

        if sacfile != None:
            self.read(sacfile)


    def read(self, sacfile):
        """
        Read sac file: append data, times, and headers
        """
        if os.path.isfile(sacfile):
            sacfile = sacfile
        else:
            sacfile = None
            print(f"Error: sac file does not exist: '{sacfile}'")
            exit(-1)

        # try read sac data and times
        try:
            st = obspy.read(sacfile, format="SAC", headonly=False)
            tr = st[0]
            data = tr.data
            times = tr.times()
        except Exception as e:
            print(e)
            print(f"Error: could not read sac data: '{sacfile}'")
            exit(-1)

        headers = {}
        for key in  tr.stats.sac.keys():
            headers[f"{key}"] = tr.stats.sac[f"{key}"]

        self.append(data, times, headers)


    def append(self, data, times, headers):

        if type(headers) != dict:
            print("Error: appended 'headers' must be a type 'dict'.")
            exit(1)

        if len(times) != len(data):
            print("Error: number of elements in 'times' and 'data' do not match!")
            exit(1)

        self.data = data
        self.times = times
        self.headers = headers
        self.update_headers()


    def update_headers(self):

        # 'knetwk', 'kstnm', 'kcmpnm', 'stla', 'stlo', 'stel'
        if 'station' in self.headers.keys():
            self.headers['kstnm'] = self.headers['station']
        if 'network' in self.headers.keys():
            self.headers['knetwk'] = self.headers['network']
        if 'channel' in self.headers.keys():
            self.headers['kcmpnm'] = self.headers['channel']
        if 'latitude' in self.headers.keys():
            self.headers['stla'] = self.headers['latitude']
        if 'longitude' in self.headers.keys():
            self.headers['stlo'] = self.headers['longitude']
        if 'elevation' in self.headers.keys():
            self.headers['stel'] = self.headers['elevation']

        for key in Required_SAC_Headers:
            if key not in self.headers.keys():
                print(f"Error: required header information is missing: '{key}'")
                exit(1)

        self.headers['knetwk'] = str(self.headers['knetwk'])
        self.headers['kstnm'] = str(self.headers['kstnm'])
        self.headers['kcmpnm'] = str(self.headers['kcmpnm'])
        self.headers['stla'] = float(self.headers['stla'])
        self.headers['stlo'] = float(self.headers['stlo'])
        self.headers['stel'] = float(self.headers['stel'])

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

        # update self.endtime
        self.endtime = self.starttime + (self.times[2] - self.times[1])

        # these two are problematic!
        if 'iftype' in self.headers.keys():
            del self.headers['iftype']
        if 'iztype' in self.headers.keys():
            del self.headers['iztype']

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

        header_keys = list(self.headers.keys())
        for key in header_keys:
            if key.upper() not in All_Possible_SAC_Headers:
                self.headers.pop(key)


    def filter(self, filter_design, domain):
        if type(filter_design) == list:
            filter_design = np.array(filter_design, dtype=float)

        if domain.lower() in ["f", "frequency", "freq"]: # apply frequency domain filter
            zero_pad = np.zeros(len(self.data) - len(filter_design)).tolist()
            filt_zero_padded = np.array(filter_design.tolist() + zero_pad)
            timeseries_filtered_fdomain = np.multiply(fft(self.data), filt_zero_padded)
            timeseries_filtered_tdomain = np.real(ifft(timeseries_filtered_fdomain))
        elif domain.lower() in ["t", "time"]: # apply time domain filter
            timeseries_filtered_tdomain = np.multiply(self.data, filter_design)
        else:
            print("Error: type choices: 'F' (Frequency Domain) or 'T' (Time Domain)")
            exit(1)
        self.data = timeseries_filtered_tdomain


    def resample(self, fs, lowpass=True):
        """
        fs: new sampling frequency
        lowpass: apply appropriate low pass filtering before signal decimation (True/False)

        """
        new_delta = 1 / fs
        new_times = np.arange(self.times[0], self.times[-1] + new_delta, new_delta)
        current_fs = 1 / self.headers['delta']
        if fs < current_fs and lowpass:
            # step 1: low pass filter
            order = 10
            nyq_freq = 0.5 * fs
            b, a = butter(order, nyq_freq / 3, btype='low', analog=False, fs=fs)
            data_lp = filtfilt(b, a, self.data)
            # step 2: linear 1D interpolation
            new_data = np.interp(new_times, self.times, data_lp)
        else:
            new_data = np.interp(new_times, self.times, self.data)
        # save changes
        self.data = new_data
        self.times = new_times
        self.update_headers()



    def demean(self):
        data_mean = np.nanmean(self.data)
        self.data = self.data - data_mean


    def detrend(self, method='spline', order=3, dspline=500):
        tr = obspy.core.trace.Trace()
        tr.data = np.array(self.data)
        tr.times = np.array(self.times)
        if method == 'spline':
            try:
                tr.detrend(method, order=order, dspline=dspline)
            except:
                tr.detrend('demean')
        elif method == 'polynomial':
            try:
                tr.detrend(method, order=order)
            except:
                tr.detrend('demean')
        else:
            tr.detrend(method)
        self.data = tr.data


    def narrowband(self, c, unit, taper_ratio=0.5, return_filter=False, plot=False):
        """
        c: central fequency/period of the narrow band filter (specified by 'unit')
        unit: period/frequency
        taper_ratio: taper_width = taper_ratio * central_frequency
        return_filter: return filter design [faxis, filter_design]
        plot: plot filter design
        """
        self.bandpass(c, c, unit=unit, taper_ratio=taper_ratio, plot=plot)


    def bandpass(self, c1, c2, unit, taper_ratio=0.1, return_filter=False, plot=False):
        """
        c1, c2: corner periods/frequencies (specified by 'unit')
        unit: period/frequency
        taper_ratio: corner_taper_width = taper_ratio * central_frequency
        plot: show filter design plot
        return_filter: return filter design [faxis, filter_design]

        **Notes**
            1) if taper_ratio==0.0, it is going to be a box bandpass filter
            1) if c1==c2, it is going to be a narrowband filter and taper_ratio must be > 0.0

        """
        # sanity checks
        if unit.lower() in ['f', 'freq', 'frequency']:
            f1 = float(c1)
            f2 = float(c2)
        elif unit.lower() in ['p', 'period']:
            f1 = float(1 / c2)
            f2 = float(1 / c1)
        else:
            print("Error: bandpass parameter: unit must be given as either 'p' or 'f'")
            exit(1)
        if f1 > f2:
            print("Error: bandpass c1 and c2: c1 must be smaller than c2 (or equal to i.e. narrowband)!")
            exit(1)
        if f1 == f2 and taper_ratio==0.0:
            print("Error: for narrowband filtering (c1=c2), 'taper_ratio' must be larger than 0.0!")
            exit(1)

        # design frequency domain bandpass filter
        fc = (f1 + f2) / 2 # central frequency
        faxis = np.arange(0, np.floor(self.headers['npts']/2) + 1, 1)\
                / self.headers['delta'] / self.headers['npts']
        # calculate filter taper width as a function of central frequency
        freq_interval = faxis[1] - faxis[0]
        taper_width = fc * taper_ratio

        if fc == f1 == f2: # narrowband filtering

            filter_design = np.exp(-1* np.power(faxis-fc, 2.) / 2. / np.power((taper_width * fc), 2))

        else: # bandpass filtering

            # first, design box filter
            filter_design = np.zeros(len(faxis))
            for ifreq, freq in enumerate(faxis):
                if freq >= f1 and freq <= f2:
                    filter_design[ifreq] = 1.0
                    max_index = ifreq
            
            # design hanning taper
            nhann = int((taper_width * 2) / freq_interval)
            if nhann % 2 == 0:
                nhann += 1
            full_hanning = np.hanning(nhann)
            half_hanning = []
            for i in range(nhann):
                half_hanning.append(full_hanning[i])
                if full_hanning[i] == 1.0:
                    break
            # left hanning taper
            ihann = 0
            nhann = int(np.floor(nhann / 2)) + 1
            for ifreq, freq in enumerate(faxis):
                if freq >= (f1 - taper_width) and ihann < nhann:
                    filter_design[ifreq] = half_hanning[ihann]
                    ihann += 1
            # right hanning taper
            for ihann in range(nhann):
                indx = max_index - ihann + int(taper_width / freq_interval)
                filter_design[indx] = half_hanning[ihann]

        if plot:
            fig = plt.figure(figsize=(10, 3))
            plt.plot(faxis, filter_design)
            plt.show()
            plt.close()

        if return_filter:
            return [faxis, filter_design]

        self.filter(filter_design, domain='frequency')



    def taper(self, width=50, tmin=None, tmax=None):
        """
        width: taper width in seconds
        tmin: minimum time that data tapers to zero
        tmax: maximum time that data tapers to zero

        """
        if tmin == None or tmin < (np.nanmin(self.times) + width):
            tmin = np.nanmin(self.times) + width
        if tmax == None or tmax > (np.nanmax(self.times) - width):
            tmax = np.nanmax(self.times) - width
        # first, design box filter
        filt = np.zeros(len(self.data))
        for itime, time in enumerate(self.times):
            if time >= tmin and time <= tmax:
                filt[itime] = 1.0
                max_index = itime
        # design hanning taper
        nhann = int((width * 2) / self.headers['delta'])
        if nhann % 2 == 0:
            nhann += 1
        full_hanning = np.hanning(nhann)
        half_hanning = []
        for i in range(nhann):
            half_hanning.append(full_hanning[i])
            if full_hanning[i] == 1.0:
                break
        # left hanning taper
        ihann = 0
        nhann = int(np.floor(nhann / 2)) + 1
        for itime, time in enumerate(self.times):
            if time >= tmin and ihann < nhann:
                filt[itime] = half_hanning[ihann]
                ihann += 1
        # right hanning taper
        for ihann in range(nhann):
            indx = max_index - ihann
            filt[indx] = half_hanning[ihann]

        self.filter(filt, domain='time')



    def plot(self, output=None):
        fig = plt.figure(figsize=(10,3))
        plt.plot(self.times, self.data)
        plt.title(f"{self.headers['knetwk']}.{self.headers['kstnm']}.{self.headers['kcmpnm']}")
        plt.xlabel('Time (s)')
        plt.ylabel('Amplitude')
        plt.tight_layout()
        # plt.ylim([-2000, 2000]) # just for test
        if output:
            fext = os.path.splitext(output)[1]
            if not fext:
                fext = 'png'
                output = f"{output}.{fext}"
            plt.savefig(output, dpi=400)
        else:
            plt.show()
        plt.close()


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


    def get_endtime(self):
        return self.endtime


