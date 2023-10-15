#!/usr/bin/env python

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.fft import fft, ifft, rfft, irfft

from core.station import Station
from core.event import Event

########FUNCTIONS########

def periods_from_frequencies(freqs):
    """
    Calculate a list of periods from a list of frequencies

    Types
    ----------
    list(float): periods, freqs

    Arguments
    ----------
    freqs: a list of frequencies in hertz

    Returns
    ----------
    periods: list of periods in seconds
    """
    periods = []
    for f in freqs:
        periods.append(1/f)
    return periods



def frequencies_from_periods(periods):
    """
    Calculate a list of periods from a list of frequencies

    Types
    ----------
    list(float): periods, freqs

    Arguments
    ----------
    periods: list of periods in seconds

    Returns
    ----------
    freqs: list of frequencies in hertz
    """
    freqs = []
    for prd in periods:
        freqs.append(1/prd)
    return freqs



def design_gaussian_filters(freqs,dt,npts,filter_width_range):
    """
    Gaussian filters for a collection of central frequencies
    to be applied in the frequency domains

    Types
    ----------
    list(float):  freqs, faxis, filter_width_range
    list(list(float)): gaussian_filters
    int: npts
    float: dt, minwidth, maxwidth

    Arguments
    ----------
    freqs: central frequencies
    dt: timeseries delta
    npts: number of points
    filter_width_range: minimum and maximum filter width range
    [minwidth, maxwidth] = filter_width_range
    minwidth: minimum filter width (fraction of central frequency)
    maxwidth: maximum filter width (fraction of central frequency)

    Returns
    ----------
    faxis: frequency axis (horizontal axis) for the designed filters
    gaussian_filters: list of designed gaussian filters (vertical axis; values in 0-1 range)
    """
    freqs = np.array(freqs)
    faxis = np.arange(0, np.floor(npts/2) + 1, 1) / dt / npts

    [minwidth, maxwidth] = filter_width_range
    filter_widths = maxwidth - (maxwidth - minwidth) /\
             np.ptp(freqs) * (freqs - np.min(freqs))
    gaussian_filters = []
    for ifreq, freq in enumerate(freqs):
        gaussian_filters.append(
            np.exp(-1* np.power(faxis-freq, 2.) / 2. / np.power((filter_widths[ifreq] * freq), 2))
        )
    return [faxis, gaussian_filters]


def apply_filter(timeseries, filt):
    """
    Applies a filter to the given timesseries in the frequency domain
    Function procedure:
        1) take timeseries to the frequency domain (fft)
        2) multiply the spectrum to the given filter
        3) back to the time domain and return filtered timeseries (ifft)

    Types
    ---------
    list(float): timeseries, filt

    Arguments
    ---------
    timeseries: timeseries in the time domain (amplitudes only)
    filt: designed filter in the frequency domain

    Returns
    ---------
    timeseries_filtered: filtered timeseries in the time domain
    """
    zero_pad = np.zeros(len(timeseries) - len(filt)).tolist()
    filt_zero_padded = np.array(filt.tolist() + zero_pad)
    timeseries_filtered_fdomain = np.multiply(fft(timeseries), filt_zero_padded)
    timeseries_filtered = ifft(timeseries_filtered_fdomain)
    return np.real(timeseries_filtered)


def calculate_narrow_band_if_good(sta_obj, filt, freq, minvel, maxvel, peak_tol):
    """
    Applies a narrow band Gaussian filter considering the 
    group velocity range and a peak tolerance value if all good!

    """
    tmin = sta_obj.headers["dist"]/maxvel + sta_obj.headers["o"];
    tmax = sta_obj.headers["dist"]/minvel + sta_obj.headers["o"];
    if tmax > np.max(sta_obj.times) or tmin < np.min(sta_obj.times):
        print(f"Error: Station '{sta_obj.headers['kstnm']}' does not contain enough data (freq={freq}; velocities:[{minvel}, {maxvel}])")
        return None
    filtered_timeseries = apply_filter(sta_obj.data, filt)
    for i, t in enumerate(sta_obj.times):
        print(t, tmin, tmax)
        if t < tmin or t > tmax:
            filtered_timeseries[i] = 0.
    return filtered_timeseries




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



##########################

#========Adjustable Parameters=========#
periods = [20, 25, 30, 40, 50, 60, 80, 100, 120, 140, 160, 180, 200]
filter_width_range = [0.06, 0.1]
events = [
    os.path.abspath(os.path.join('test_events','08174235631')),
    os.path.abspath(os.path.join('test_events','08178211916'))
]
maxgroupv = 2.5
mingroupv = 5


test_event = events[1]
test_sacfile = os.path.join(events[1], "08178211916_117A.LHZ")
#======================================#


if __name__=="__main__":


    # periods/frequencies
    periods = sorted(periods)
    freqs = frequencies_from_periods(periods)
    nfreqs = len(freqs)

    # test station read
    #------------------
    sta = Station(test_sacfile)
    fax, gaussian_filters = design_gaussian_filters(freqs,
                            sta.headers['delta'],
                            len(sta.data),
                            filter_width_range)


    # plot designed gaussian filters
    # ------------------------------
    # fig = plt.figure(figsize=(10,5))
    # ax = fig.add_subplot(111)
    # for i in range(nfreqs):
    #     plt.plot(fax, gaussian_filters[i])
    # ax.set_title("Gaussian Filters")
    # plt.xlim([0, 0.3])

    # plot original versus narrow-band Gaussian filtered
    # --------------------------------------------------
    # fig, ax = plt.subplots(nfreqs+1)
    # fig.set_figwidth(6*0.8)
    # fig.set_figheight(8*0.8)
    # ax[0].plot(sta.times, sta.data, label='Unfiltered')
    # ax[0].legend(loc='upper right')
    # ax[0].set_xticks([],[])
    # ax[0].set_yticks([],[])
    # for i in range(nfreqs):
    #     filt = gaussian_filters[i]
    #     filtered = apply_filter(sta.data, filt)
    #     ax[i+1].plot(sta.times, filtered, label=f'{periods[i]} s')
    #     ax[i+1].legend(loc='upper right')
    #     ax[i+1].set_xticks([],[])
    #     ax[i+1].set_yticks([],[])
    # plt.tight_layout()

    # plot original versus final filtered seismograms
    # -----------------------------------------------
    tmin = sta.headers["dist"]/maxgroupv + sta.headers["o"];
    tmax = sta.headers["dist"]/mingroupv + sta.headers["o"];
    if tmax > np.max(sta.times) or tmin < np.min(sta.times):
        print(f"Error: Station '{sta.headers['kstnm']}' does not contain enough data")
        exit(-1)

    fig, ax = plt.subplots(nfreqs+1)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(8*0.8)
    ax[0].plot(sta.times, sta.data, label='Unfiltered')
    ax[0].scatter([tmin, tmax],[0, 0], c='red')
    ax[0].legend(loc='upper right')
    ax[0].set_xticks([],[])
    ax[0].set_yticks([],[])
    for i in range(nfreqs):
        filt = gaussian_filters[i]
        filtered = calculate_narrow_band_if_good(sta, filt, freqs[i], mingroupv, maxgroupv, 0.1)
        ax[i+1].plot(sta.times, filtered, label=f'{periods[i]} s')
        ax[i+1].legend(loc='upper right')
        ax[i+1].set_xticks([],[])
        ax[i+1].set_yticks([],[])
    plt.tight_layout()

    plt.show()
    exit(0)




    # test events
    evt = Event(extensions=["LHZ"])
    evt.read_folder(test_event)
    # evt.read_folder(test_event2, extensions=["LHZ"])


    # design gaussian filters
    fax, gaussian_filters = design_gaussian_filters(freqs, 1, 3000, 0.06, 0.10)




