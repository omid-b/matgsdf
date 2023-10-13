#!/usr/bin/env python

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from core.station import Station
from core.event import Event

########FUNCTIONS########

def periods_from_frequencies(freqs):
    """
    Calculate a sorted list of periods from a list of frequencies

    Types
    ----------
    list(float): periods, freqs

    Arguments
    ----------
    freqs: a list of frequencies in hertz

    Returns
    ----------
    periods: a sorted list of periods in seconds
    """
    periods = []
    for f in freqs:
        periods.append(1/f)
    return sorted(periods)



def frequencies_from_periods(periods):
    """
    Calculate a sorted list of periods from a list of frequencies

    Types
    ----------
    list(float): periods, freqs

    Arguments
    ----------
    periods: list of periods in seconds

    Returns
    ----------
    freqs: sorted list of frequencies in hertz
    """
    freqs = []
    for prd in periods:
        freqs.append(1/prd)
    return sorted(freqs)



def design_gaussian_filters(freqs,dt,npts,minwidth,maxwidth):
    """
    Gaussian filters for a collection of central frequencies
    to be applied in the frequency domains

    Types
    ----------
    list(float):  freqs, faxis
    list(list(float)): gaussian_filters
    int: npts
    float: dt, minwidth, maxwidth

    Arguments
    ----------
    freqs: central frequencies
    dt: timeseries delta
    npts: number of points
    minwidth: minimum filter width (fraction of central frequency)
    maxwidth: maximum filter width (fraction of central frequency)

    Returns
    ----------
    faxis: frequency axis (horizontal axis) for the designed filters
    gaussian_filters: list of designed gaussian filters (vertical axis; values in 0-1 range)
    """
    freqs = np.array(freqs)
    faxis = np.arange(0, np.floor(npts/2) + 1, 1) / dt / npts
    widths = maxwidth - (maxwidth - minwidth) /\
             np.ptp(freqs) * (freqs - np.min(freqs))
    gaussian_filters = []
    for ifreq, freq in enumerate(freqs):
        gaussian_filters.append(
            np.exp(-1* np.power(faxis-freq, 2.) / 2. / np.power((widths[ifreq] * freq), 2))
        )
    return [faxis, gaussian_filters]


def apply_filter(timeseries, filter):
    """
    Applies a filter to the given timesseries in the frequency domain
    Function procedure:
        1) take timeseries to the frequency domain (fft)
        2) multiply the spectrum to the given filter
        3) back to the time domain and return filtered timeseries (ifft)

    Types
    ---------
    list(float): timeseries, filter

    Arguments
    ---------
    timeseries: timeseries in the time domain (amplitudes only)
    filter: designed filter in the frequency domain (values between 0 and 1)

    Returns
    ---------
    timeseries_filt: filtered timeseries in the time domain
    """
    pass



def groupv_fit(dist,time,snr,mingroupv,maxgroupv):
    """
    Apply weighted least square polynomial fitting to the data

    Types
    -------
    XXX

    Arguments
    -------
    XXX

    Returns
    -------
    XXX
    
    """
    pass



def read_event_sacfiles(event_dir, file_extensions=["sac", "SAC"], verbose=True):
    """
    Given and event directory containing sac files, this function
    first check the sac file headers and then returns an event object
    that is a dictionary that contains all the necessary informations
    for that event.

    Types
    ------------
    str: event_dir
    list(str): file_extensions
    bool: verbose

    Arguments
    ------------
    event_dir: full path to directory contining the sac files
    file_extensions: list of sac file extensions
    verbose: print all the errors and warnings

    Returns
    ------------
    event: an event object i.e. a python dictionary:
        event = {
            "evla":  <float: event_latitude>,  # event latitude
            "evlo":  <float: event_longitude>, # event longitude
            "evlo":  <float: event_longitude>, # event longitude
            "otime":     <float: event_origin_time>,   # event origin time (float value in time axis)
            "otime_str": <str: event_origin_time>,     # event origin time in string in YYJJJHHMMSS format
            "otime_utc": <obspy UTCDateTime: event_origin_time>, # event origin time obspy UTCDateTime object
            "times":     <list(float): times>, # time axis of the timeseries
            "data":      <list(list(float)): station_data>, # list of time series data
        }

    """
    pass

##########################

#========Adjustable Parameters=========#
periods = [20, 25, 32, 40, 50, 60, 80, 100]
events = [
    os.path.abspath(os.path.join('test_events','08174235631')),
    os.path.abspath(os.path.join('test_events','08178211916'))
]
#======================================#


if __name__=="__main__":

    test_sacfile = os.path.join("test_events", "08174235631", "08174235631_113A.LHZ")
    test_event = os.path.join("test_events", "08174235631")

    # test station
    # sta = Station(test_sacfile)
    # print(sta.get_headers())

    # test events
    evt = Event()

    evt.read_folder(test_event, extensions=["LHZ"])



    exit(0)

    # periods/frequencies
    periods = sorted(periods)
    freqs = frequencies_from_periods(periods)
    nfreqs = len(freqs)

    # design gaussian filters
    fax, gaussian_filters = design_gaussian_filters(freqs, 1, 3000, 0.06, 0.10)

    # plot designed gaussian filters
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    for i in range(nfreqs):
        plt.plot(fax, gaussian_filters[i])
    ax.set_title("Gaussian Filters")
    plt.xlim([0, 0.3])
    plt.show()



