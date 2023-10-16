#!/usr/bin/env python

import os
import sys
import math


from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from core.event import Event
from core.station import Station
from core.xcorr import XCorr

from utils.calc import frequencies_from_periods
from utils.calc import periods_from_frequencies

from filtering.isolation import apply_filter
from filtering.isolation import design_gaussian_filters
from filtering.isolation import calculate_narrow_band_if_good


#========Adjustable Parameters=========#
periods = [20, 25, 30, 40, 50, 60, 80, 100, 120, 140, 160]
filter_width_range = [0.06, 0.1]
events = [
    os.path.abspath(os.path.join('_datasets', 'test_events','08174235631')),
    os.path.abspath(os.path.join('_datasets', 'test_events','08178211916'))
]
mingroupv = 2.5
maxgroupv = 5

test_egfs = os.path.join("_datasets", "test_egfs")
test_egfs2 = os.path.join("_datasets", "test_egfs_2")
test_event = events[1]
test_sacfile = os.path.join(events[1], "08178211916_117A.LHZ")

# periods/frequencies
periods = sorted(periods)
freqs = frequencies_from_periods(periods)
nfreqs = len(freqs)

# test station read
#------------------
sta = Station(test_sacfile)
#======================================#

###### TEST FUNCTIONS #######

def test_plot_egfs():
    # read and plot a EGFs (XCorr object)
    #-----------------------------
    egfs = XCorr(test_egfs)
    print("----test_egfs----")
    print("cmp: ", egfs.component)
    print("tmin: ", egfs.tmin)
    print("tmax: ", egfs.tmax)
    print("delta: ", egfs.delta)
    print("min_dist: ", egfs.min_dist)
    print("max_dist: ", egfs.max_dist)
    print("avg_dist: ", egfs.avg_dist)
    print("sym: ", egfs.sym)

    nsta = len(egfs.stations)
    fig, ax = plt.subplots(nsta)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(10*0.8)
    for ista in range(nsta):
        ax[ista].plot(egfs.stations[ista].times,\
                   egfs.stations[ista].data,\
                   label="%.2f km" %(egfs.stations[ista].headers['dist']))
        ax[ista].legend(loc="upper right")
        ax[ista].set_yticks([],[])
        ax[ista].set_xlim([-500, 500]) # test_egfs
    ax[0].set_title("EGFs sorted by distance")

    egfs = XCorr(test_egfs2, component="ZZ")
    print("----test_egfs2----")
    print("cmp: ", egfs.component)
    print("tmin: ", egfs.tmin)
    print("tmax: ", egfs.tmax)
    print("delta: ", egfs.delta)
    print("min_dist: ", egfs.min_dist)
    print("max_dist: ", egfs.max_dist)
    print("avg_dist: ", egfs.avg_dist)
    print("sym: ", egfs.sym)
    nsta = len(egfs.stations)
    fig, ax = plt.subplots(nsta)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(10*0.8)
    for ista in range(nsta):
        ax[ista].plot(egfs.stations[ista].times,\
                   egfs.stations[ista].data,\
                   label="%.2f km" %(egfs.stations[ista].headers['dist']))
        ax[ista].legend(loc="upper right")
        ax[ista].set_yticks([],[])
        # ax[ista].set_xlim([-500, 500]) # test_egfs
        ax[ista].set_xlim([0, 5]) # test_egfs2
    ax[0].set_title("EGFs sorted by distance")
    plt.show()


def test_plot_event():
    # read and plot a single event
    #-----------------------------
    event = Event(test_event, extensions=["LHZ"], sort_by_dist=True)
    # nsta = len(event.stations)
    nsta = len(event.stations) - 25 # just for test!
    fig, ax = plt.subplots(nsta)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(10*0.8)
    for ista in range(nsta):
        ax[ista].plot(event.stations[ista].times,\
                   event.stations[ista].data,\
                   label="%.2f km" %(event.stations[ista].headers['dist']))
        ax[ista].legend(loc="upper right")
        ax[ista].set_yticks([],[])
        # ax[ista].set_xlim([500, 1000])
    ax[0].set_title("Event data sorted by distance")
    plt.show()


    
    
def test_narrow_band_gaussians():
    # design and plot narrow band gaussian filters
    # ------------------------------
    fax, gaussian_filters = design_gaussian_filters(freqs,
                            sta.headers['delta'],
                            len(sta.data),
                            filter_width_range)
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    for i in range(nfreqs):
        plt.plot(fax, gaussian_filters[i], label=f"{periods[i]} s")
    plt.legend(loc="best")  
    plt.xlabel("Frequency (hertz)")  
    plt.ylabel("Filter Amplitude")  
    ax.set_title("Gaussian Filters")
    plt.xlim([0, 0.3])
    plt.tight_layout()
    plt.show()


def test_narrow_band_filtered():
    fax, gaussian_filters = design_gaussian_filters(freqs,
                            sta.headers['delta'],
                            len(sta.data),
                            filter_width_range)
    # plot original versus narrow-band Gaussian filtered
    # --------------------------------------------------
    fig, ax = plt.subplots(nfreqs+1)
    fig.set_figwidth(6*0.8)
    fig.set_figheight(8*0.8)
    ax[0].plot(sta.times, sta.data, label='Unfiltered')
    ax[0].legend(loc='upper right')
    ax[0].set_xticks([],[])
    ax[0].set_yticks([],[])
    for i in range(nfreqs):
        filt = gaussian_filters[i]
        filtered = apply_filter(sta.data, filt)
        ax[i+1].plot(sta.times, filtered, label=f'{periods[i]} s')
        ax[i+1].legend(loc='upper right')
        ax[i+1].set_xticks([],[])
        ax[i+1].set_yticks([],[])
    plt.tight_layout()
    plt.show()


def test_isolated_filtered():
    fax, gaussian_filters = design_gaussian_filters(freqs,
                            sta.headers['delta'],
                            len(sta.data),
                            filter_width_range)
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
        filtered = calculate_narrow_band_if_good(sta, filt, freqs[i], mingroupv, maxgroupv, 1)
        ax[i+1].plot(sta.times, filtered, label=f'{periods[i]} s')
        ax[i+1].legend(loc='upper right')
        ax[i+1].set_xticks([],[])
        ax[i+1].set_yticks([],[])
    plt.tight_layout()
    plt.show()



if __name__=="__main__":
    test_plot_egfs()
    test_plot_event()
    test_narrow_band_gaussians()
    test_narrow_band_filtered()
    test_isolated_filtered()




