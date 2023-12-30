#!/usr/bin/env python

# -*- coding: utf-8 -*-
r"""
Seismic SAC file trace class
==============================

read, process, and output processed single seismic traces

Class Trace: single station SAC format trace

Coded by: omid.bagherpur@gmail.com; github: omid-b
Update: 2023-12-30 (WORK in Progress!)

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

from scipy.fft import fft
from scipy.fft import ifft
from scipy.fft import rfft
from scipy.fft import irfft

from scipy.signal import butter
from scipy.signal import freqz

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
        self.response_removed = False

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
            raise ValueError(f"SAC file does not exist: '{sacfile}'")

        # try read sac data and times
        try:
            st = obspy.read(sacfile, format="SAC", headonly=False)
            tr = st[0]
            data = tr.data
            times = tr.times()
        except Exception:
            raise ValueError(f"Could not read sac data: '{sacfile}'")

        headers = {}
        for key in  tr.stats.sac.keys():
            headers[f"{key}"] = tr.stats.sac[f"{key}"]

        self.append(data, times, headers)


    def append(self, data, times, headers):
        """
        Set the trace data, times and headers by appending the lists and headers dictionary

        Note: this is an alternative to self.read() method

        """

        if type(headers) != dict:
            raise ValueError("Appended 'headers' must be a type 'dict'.")

        if len(times) != len(data):
            raise ValueError("Number of elements in 'times' and 'data' do not match!")

        self.data = data
        self.times = times
        self.headers = headers
        self.update_headers()


    def update_headers(self):
        """
        Class method to update trace headers after processing

        """

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
                raise ValueError(f"Required header information is missing: '{key}'")

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



    def apply_filter(self, filter_design, domain):
        """
        Apply time or frequency domain filter to the trace data

        filter_design: filter design, values between 0 and 1
        domain: frequency or time; choices: ["f", "frequency", "freq", "t", "time"]

        """
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
            raise ValueError("Type choices: 'F' (Frequency Domain) or 'T' (Time Domain)")

        self.data = timeseries_filtered_tdomain



    def demean(self):
        """
        Demean trace: subtract the mean from the data

        """
        data_mean = np.nanmean(self.data)
        self.data = self.data - data_mean


    def detrend(self, method='spline', order=3, dspline=500):
        """
        Detrend trace using various methods
        
        method: detrending method; choices: ['spline', 'polynomial', 'linear', 'demean']
        order: spline or polynomial order
        dspline: number of dsplines (spline method only)

        """
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


    def taper(self, width=50, tmin=None, tmax=None):
        """
        Apply taper to the tace in time range: [tmin, tmax]

        width: taper width in seconds
        tmin: minimum time that data tapers to zero
        tmax: maximum time that data tapers to zero

        **Note**
        If 'tmin' and 'tmax' not given, they will be set automatically
        to the begin and end of the trace (similar to SAC command)

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

        self.apply_filter(filt, domain='time')


    def cut(self, cut_begin, cut_end):
        """
        cut trace in range: [cut_begin, cut_end]

        cut_begin: cut begin time
        cut_end: cut end time

        """
        current_fs = 1 / self.headers['delta']
        new_times = np.arange(cut_begin, cut_end + current_fs, current_fs)
        new_data = np.zeros(len(new_times))
        new_times_interp = []
        for t in new_times:
            if t >= self.headers['b'] and t <= self.headers['e']:
                new_times_interp.append(t)
        # interpolate at times in new_times_interp (the rest is set to zero; zerofill)
        new_data_interp = np.interp(new_times_interp, self.times, self.data)
        for i, t in enumerate(new_times):
            if t in new_times_interp:
                indx = new_times_interp.index(t)
                new_data[i] = new_data_interp[indx]
        self.times = new_times
        self.data = new_data
        self.update_headers()


    def whiten(self):
        pass
        # XXX


    def normalize(self, type='onebit'):
        pass
        # XXX


    def calc_snr(self, signal_window=[], noise_window=[], velocities=[]):
        pass
        # XXX


    def remove_response(self, xml_file, output='VEL', water_level=60, 
                        pre_filt=None, update_headers=True):
        try:
            inv = obspy.read_inventory(xml_file)
            xml_knetwk = inv[0].code.split()[0]
            xml_kstnm  = inv[0][0].code.split()[0]
            xml_kcmpnm = inv[0][0][0].code.split()[0]
            xml_stla   = float(inv[0][0].latitude)
            xml_stlo   = float(inv[0][0].longitude)
            xml_stel   = float(inv[0][0].elevation)
            xml_cmpaz  = float(inv[0][0][0].azimuth)
            xml_cmpinc = float(inv[0][0][0].dip) + 90
        except:
            raise ValueError("Could not read input XML file.")

        if  xml_knetwk != self.headers['knetwk'] or \
            xml_kstnm != self.headers['kstnm'] or \
            xml_kcmpnm != self.headers['kcmpnm']:
            raise ValueError('XML info tag (i.e. "network.station.component") does not match')
        try:
            tr = SACTrace( data=np.array(self.data), **self.headers)
            tr.remove_response(inventory=inv, output=output, water_level=water_level, pre_filt=pre_filt)
        except:
            raise ValueError('Could not remove instrument response from data')
        self.response_removed = True


    def resample(self, fs, lowpass=True):
        """
        Resample trace at a different sampling frequency

        fs: new sampling frequency
        lowpass: apply appropriate low pass filtering before signal decimation (True/False)

        """
        new_delta = 1 / fs
        new_times = np.arange(self.times[0], self.times[-1] + new_delta, new_delta)
        current_fs = 1 / self.headers['delta']
        if fs < current_fs and lowpass:
            # step 1: low pass filter
            nyq_freq = 0.5 * fs
            self.lowpass(nyq_freq / 2, 'f', type='butter', butter_order=10)

        # step 2: linear 1D interpolation
        new_data = np.interp(new_times, self.times, self.data)
        # save changes
        self.data = new_data
        self.times = new_times
        self.update_headers()



    def lowpass(self, c, unit, type='butter', butter_order=3, hann_width=None, plot=False, return_filter=False):
        """
        Apply lowpass filter to the trace

        c: cutoff period/frequency (specified by 'unit')
        unit: period/frequency; ; choices: ['f', 'freq','frequency', 'p', 'period']
        type: Hanning or Butterworth taper?; choices: ['hann', 'butter']
        butter_order: Butterworth filter tapering order
        hann_width (only if type='hann''): if hann_width=None, hann_width = 0.5 * cutoff_frequency
        plot: show filter design plot
        return_filter: return frequency axis & filter design arrays [faxis, filter_design]

        """
        # sanity checks
        if type.lower() not in ['hann', 'butter']:
            raise ValueError("Type must be either 'butter' or 'hann'!")

        if unit.lower() in ['f', 'freq', 'frequency']:
            f = float(c)
        elif unit.lower() in ['p', 'period']:
            f = float(1 / c)
        else:
            raise ValueError("Lowpass filter parameter: unit must be either 'p' or 'f'")

        # design and apply filter

        faxis = np.arange(0, np.floor(self.headers['npts']/2) + 1, 1)\
                    / self.headers['delta'] / self.headers['npts']

        if type == 'butter':

            current_fs = 1 / self.headers['delta']
            b, a = butter(butter_order, f, btype='lowpass', analog=False, fs=current_fs)
            _, filter_design = freqz(b, a, fs=current_fs, worN=len(faxis))
            filter_design = np.absolute(filter_design)
            
        elif type == 'hann':

            # first, design box filter
            filter_design = np.zeros(len(faxis))
            for ifreq, freq in enumerate(faxis):
                if freq <= f:
                    filter_design[ifreq] = 1.0
                    max_index = ifreq
            
            # design lowpass hanning taper
            if hann_width == None:
                hann_width = 0.5 * f
            freq_interval = faxis[1] - faxis[0]
            nhann = int((hann_width * 2) / freq_interval)
            if nhann % 2 == 0:
                nhann += 1
            full_hanning = np.hanning(nhann)
            half_hanning = []
            for i in range(nhann):
                half_hanning.append(full_hanning[i])
                if full_hanning[i] == 1.0:
                    break
            ihann = 0
            nhann = int(np.floor(nhann / 2)) + 1
            # lowpass filter => right hanning taper
            for ihann in range(nhann):
                indx = max_index - ihann + int(hann_width / freq_interval)
                filter_design[indx] = half_hanning[ihann]

        self.apply_filter(filter_design, 'frequency')

        if plot:
            fig = plt.figure(figsize=(10, 3))
            plt.plot(faxis, filter_design)
            plt.show()
            plt.close()

        if return_filter:
            return [faxis, filter_design]



    def highpass(self, c, unit, type='butter', butter_order=3, hann_width=None, plot=False, return_filter=False):
        """
        Apply highpass filter to the trace

        c: cutoff period/frequency (specified by 'unit')
        unit: period/frequency; ; choices: ['f', 'freq','frequency', 'p', 'period']
        type: Hanning or Butterworth taper?; choices: ['hann', 'butter']
        butter_order: Butterworth filter tapering order
        hann_width (only if type='hann''): if hann_width=None, hann_width = 0.5 * cutoff_frequency
        plot: show filter design plot
        return_filter: return frequency axis & filter design arrays [faxis, filter_design]

        """
        # sanity checks
        if type.lower() not in ['hann', 'butter']:
            raise ValueError("Type must be either 'butter' or 'hann'!")

        if unit.lower() in ['f', 'freq', 'frequency']:
            f = float(c)
        elif unit.lower() in ['p', 'period']:
            f = float(1 / c)
        else:
            raise ValueError("Highpass filter parameter: unit must be either 'p' or 'f'")

        # design and apply filter

        faxis = np.arange(0, np.floor(self.headers['npts']/2) + 1, 1)\
                    / self.headers['delta'] / self.headers['npts']

        if type == 'butter':

            current_fs = 1 / self.headers['delta']
            b, a = butter(butter_order, f, btype='highpass', analog=False, fs=current_fs)
            _, filter_design = freqz(b, a, fs=current_fs, worN=len(faxis))
            filter_design = np.absolute(filter_design)

        elif type == 'hann':

            # first, design box filter
            filter_design = np.zeros(len(faxis))
            for ifreq, freq in enumerate(faxis):
                if freq >= f:
                    filter_design[ifreq] = 1.0
                    max_index = ifreq
            
            # design lowpass hanning taper
            if hann_width == None:
                hann_width = 0.5 * f
            freq_interval = faxis[1] - faxis[0]
            nhann = int((hann_width * 2) / freq_interval)
            if nhann % 2 == 0:
                nhann += 1
            full_hanning = np.hanning(nhann)
            half_hanning = []
            for i in range(nhann):
                half_hanning.append(full_hanning[i])
                if full_hanning[i] == 1.0:
                    break
            ihann = 0
            nhann = int(np.floor(nhann / 2)) + 1
            # highpass filter => left hanning taper
            for ifreq, freq in enumerate(faxis):
                if freq >= (f - hann_width) and ihann < nhann:
                    filter_design[ifreq] = half_hanning[ihann]
                    ihann += 1

        self.apply_filter(filter_design, 'frequency')    

        if plot:
            fig = plt.figure(figsize=(10, 3))
            plt.plot(faxis, filter_design)
            plt.show()
            plt.close()

        if return_filter:
            return [faxis, filter_design]




    def bandstop(self, c1, c2, unit, type='butter', butter_order=3, hann_width=None, plot=False, return_filter=False):
        """
        Apply bandstop filter to the trace

        c1, c2: corner/cutoff periods/frequencies (specified by 'unit')
        unit: period/frequency; ; choices: ['f', 'freq','frequency', 'p', 'period']
        type: Hanning or Butterworth taper?; choices: ['hann', 'butter']
        butter_order: Butterworth filter tapering order
        hann_width (only if type='hann''): if hann_width=None, hann_width = 0.5 * central_frequency
        plot: show filter design plot
        return_filter: return frequency axis & filter design arrays [faxis, filter_design]

        **Notes**
            1) if hann_width==0.0, it is going to be a box bandpass filter (band-stop filter)
            1) if c1==c2, it is going to be a narrowband filter and hann_width must be > 0.0

        """
        # sanity checks
        if type.lower() not in ['hann', 'butter']:
            raise ValueError("Type must be either 'butter' or 'hann'!")

        if unit.lower() in ['f', 'freq', 'frequency']:
            f1 = float(c1)
            f2 = float(c2)
        elif unit.lower() in ['p', 'period']:
            f1 = float(1 / c2)
            f2 = float(1 / c1)
        else:
            raise ValueError("Bandstop parameter: unit must be either 'p' or 'f'")

        if f1 >= f2:
            raise ValueError("Bandstop c1 and c2: c1 must be smaller than c2!")


        # Design and Apply Filter
        faxis = np.arange(0, np.floor(self.headers['npts']/2) + 1, 1)\
                    / self.headers['delta'] / self.headers['npts']
        # calculate hanning taper width as a function of central frequency
        fc = (f1 + f2) / 2 # central frequency
        if hann_width == None:
            hann_width = 0.5 * fc

        if type=='butter':

            current_fs = 1 / self.headers['delta']
            b, a = butter(butter_order, [f1, f2], btype='bandstop', analog=False, fs=current_fs)
            _, filter_design = freqz(b, a, fs=current_fs, worN=len(faxis))
            filter_design = np.absolute(filter_design)

        elif type=='hann':

            _, fdesign_lp = self.lowpass(f1, 'f', type='hann', hann_width=hann_width, plot=False, return_filter=True)
            _, fdesign_hp = self.highpass(f2, 'f', type='hann', hann_width=hann_width, plot=False, return_filter=True)
            filter_design = np.zeros(len(faxis))
            for i, _ in enumerate(fdesign_lp):
            
                if fdesign_lp[i] > 0 and fdesign_hp[i] > 0:
                    filter_design[i] = np.max([fdesign_lp[i], fdesign_hp[i]])
                elif fdesign_lp[i] > 0 and fdesign_hp[i] == 0:
                    filter_design[i] = fdesign_lp[i]
                elif fdesign_hp[i] > 0 and fdesign_lp[i] == 0:
                    filter_design[i] = fdesign_hp[i]

        self.apply_filter(filter_design, 'frequency')    
            

        if plot:
            fig = plt.figure(figsize=(10, 3))
            plt.plot(faxis, filter_design)
            plt.show()
            plt.close()

        if return_filter:
            return [faxis, filter_design]


    def bandpass(self, c1, c2, unit, type='butter', butter_order=3, hann_width=None, plot=False, return_filter=False):
        """
        Apply bandpass filter to the trace

        c1, c2: corner period/frequency (specified by 'unit')
        unit: period/frequency; ; choices: ['f', 'freq','frequency', 'p', 'period']
        type: Hanning or Butterworth taper?; choices: ['hann', 'butter']
        butter_order: Butterworth filter tapering order
        hann_width (only if type='hann''): if hann_width=None, hann_width = 0.5 * central_frequency
        plot: show filter design plot
        return_filter: return frequency axis & filter design arrays [faxis, filter_design]

        **Notes**
            1) if hann_width==0.0, it is going to be a box bandpass filter (band-stop filter)
            1) if c1==c2, it is going to be a narrowband filter and hann_width must be > 0.0

        """
        # sanity checks
        if type.lower() not in ['hann', 'butter']:
            raise ValueError("Type must be either 'butter' or 'hann'!")

        if unit.lower() in ['f', 'freq', 'frequency']:
            f1 = float(c1)
            f2 = float(c2)
        elif unit.lower() in ['p', 'period']:
            f1 = float(1 / c2)
            f2 = float(1 / c1)
        else:
            raise ValueError("Bandpass parameter: unit must be either 'p' or 'f'")

        if f1 > f2:
            raise ValueError("Bandpass c1 and c2: c1 must be smaller than c2 (or equal to i.e. narrowband)!")

        if f1 == f2 and type=='hann' and hann_width==0.0:
            raise ValueError("For narrowband filtering (c1=c2), 'hann_width' must be larger than 0.0!")

        # Design and Apply Filter
        faxis = np.arange(0, np.floor(self.headers['npts']/2) + 1, 1)\
                    / self.headers['delta'] / self.headers['npts']
        # calculate hanning taper width as a function of central frequency
        fc = (f1 + f2) / 2 # central frequency
        if hann_width == None:
            hann_width = 0.5 * fc

        if type=='butter':
            current_fs = 1 / self.headers['delta']
            b, a = butter(butter_order, [f1, f2], btype='bandpass', analog=False, fs=current_fs)
            _, filter_design = freqz(b, a, fs=current_fs, worN=len(faxis))
            filter_design = np.absolute(filter_design)
            self.apply_filter(filter_design, 'frequency')    
        elif type=='hann':
            _, fdesign_hp = self.highpass(f1, 'f', type='hann', hann_width=hann_width, plot=False, return_filter=True)
            _, fdesign_lp = self.lowpass(f2, 'f', type='hann', hann_width=hann_width, plot=False, return_filter=True)
            filter_design = np.multiply(fdesign_hp, fdesign_lp)

        if plot:
            fig = plt.figure(figsize=(10, 3))
            plt.plot(faxis, filter_design)
            plt.show()
            plt.close()

        if return_filter:
            return [faxis, filter_design]



    def narrowband(self, c, unit, hann_width=None, plot=False, return_filter=False):
        """
        Apply narrowband filter to the trace

        c: central period/frequency (specified by 'unit')
        unit: period/frequency; ; choices: ['f', 'freq','frequency', 'p', 'period']
        hann_width: if hann_width=None, hann_width = 0.5 * central_frequency
        plot: show filter design plot
        return_filter: return frequency axis & filter design arrays [faxis, filter_design]

        **Notes**

        This narrowband filter is a bandpass type='hann' filter with c1==c2

        """
        faxis, filter_design = self.bandpass(c, c, unit=unit, type='hann', hann_width=hann_width,\
                                             plot=plot, return_filter=True)
        if return_filter:
            return return_filter



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


