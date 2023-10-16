#!/usr/bin/env python3

"""
Module for narrow-band Gaussian isolation filtering

Coded by: omid.bagherpur@gmail.com
Update: 16 Oct 2023

"""

import numpy as np
from scipy.fft import fft, ifft, rfft, irfft
# from scipy.signal import hilbert
# from scipy.signal import find_peaks


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

    envelop = filtered_timeseries.copy()
    envelop_indx = []
    for i, t in enumerate(sta_obj.times):
        if t < tmin or t > tmax:
            envelop[i] = 0.
        else:
            envelop_indx.append(i)
    envelop = envelop / np.max(envelop)
    
    envelop_diff = np.diff(envelop)

    peaks_indx = []
    peaks_vals = []
    for i in range(1, len(envelop_indx)):
        if  envelop_diff[envelop_indx[i]-1] > 0 \
            and envelop_diff[envelop_indx[i]] < 0:
            peaks_indx.append(envelop_indx[i])
            peaks_vals.append(envelop[envelop_indx[i]])

    max_peak_indx = peaks_indx[np.where(peaks_vals==np.max(peaks_vals))[0][0]]
    max_peak_val = peaks_vals[np.where(peaks_vals==np.max(peaks_vals))[0][0]]

    peaks_indx_v2 = []
    peaks_vals_v2 = []
    peak_tol = 0.01
    for i in range(1, len(peaks_indx)):
        if  envelop[peaks_indx[i]] > max_peak_val * peak_tol:
            peaks_indx_v2.append(envelop_indx[i])
            peaks_vals_v2.append(envelop[envelop_indx[i]])
    

    test = envelop.copy()
    for i in peaks_indx_v2:
        test[i] = -1
    return test


def calculate_isolation_filters(sta_obj, designed_filters, frequencies,\
                                min_group_vel, max_group_vel):
    """
    Calculates isolated timeseries
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


    