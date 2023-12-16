#!/usr/bin/env python

import os
import sys
import math


from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

from core.event import Event
from core.trace import Sac
from core.xcorr import XCorr

from utils.calc import frequencies_from_periods
from utils.calc import periods_from_frequencies

from filtering.isolation import apply_filter_fdomain
from filtering.isolation import apply_filter_tdomain
from filtering.isolation import calculate_narrow_band_if_good
from filtering.isolation import design_fdomain_gaussians
from filtering.isolation import design_tdomain_taper

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
# sac = Sac(test_sacfile)

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

    nsta = len(egfs.sacs)
    fig, ax = plt.subplots(nsta)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(10*0.8)
    for ista in range(nsta):
        ax[ista].plot(egfs.sacs[ista].times,\
                   egfs.sacs[ista].data,\
                   label="%.2f km" %(egfs.sacs[ista].headers['dist']))
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
    nsta = len(egfs.sacs)
    fig, ax = plt.subplots(nsta)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(10*0.8)
    for ista in range(nsta):
        ax[ista].plot(egfs.sacs[ista].times,\
                   egfs.sacs[ista].data,\
                   label="%.2f km" %(egfs.sacs[ista].headers['dist']))
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
    # nsta = len(event.sacs)
    nsta = len(event.sacs) - 25 # just for test!
    fig, ax = plt.subplots(nsta)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(10*0.8)
    for ista in range(nsta):
        ax[ista].plot(event.sacs[ista].times,\
                   event.sacs[ista].data,\
                   label="%.2f km" %(event.sacs[ista].headers['dist']))
        ax[ista].legend(loc="upper right")
        ax[ista].set_yticks([],[])
        # ax[ista].set_xlim([500, 1000])
    ax[0].set_title("Event data sorted by distance")
    plt.show()


    
    
def test_narrow_band_gaussian_filter_design():
    # design and plot narrow band gaussian filters
    # ------------------------------
    fax, gaussian_filters = design_fdomain_gaussians(freqs,
                            sac.headers['delta'],
                            len(sac.data),
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


def test_narrow_band_filtered_seisomgrams():
    fax, gaussian_filters = design_fdomain_gaussians(freqs,
                            sac.headers['delta'],
                            len(sac.data),
                            filter_width_range)
    # plot original versus narrow-band Gaussian filtered
    # --------------------------------------------------
    fig, ax = plt.subplots(nfreqs+1)
    fig.set_figwidth(6*0.8)
    fig.set_figheight(8*0.8)
    ax[0].plot(sac.times, sac.data, label='Unfiltered')
    ax[0].legend(loc='upper right')
    ax[0].set_xticks([],[])
    ax[0].set_yticks([],[])
    for i in range(nfreqs):
        filt = gaussian_filters[i]
        filtered = apply_filter_fdomain(sac.data, filt)
        ax[i+1].plot(sac.times, filtered, label=f'{periods[i]} s')
        ax[i+1].legend(loc='upper right')
        ax[i+1].set_xticks([],[])
        ax[i+1].set_yticks([],[])
    plt.tight_layout()
    plt.show()


# def test_isolated_filtered():
#     fax, gaussian_filters = design_fdomain_gaussians(freqs,
#                             sac.headers['delta'],
#                             len(sac.data),
#                             filter_width_range)
#     # plot original versus final filtered seismograms
#     # -----------------------------------------------
#     tmin = sac.headers["dist"]/maxgroupv + sac.headers["o"];
#     tmax = sac.headers["dist"]/mingroupv + sac.headers["o"];
#     if tmax > np.max(sac.times) or tmin < np.min(sac.times):
#         print(f"Error: Station '{sac.headers['kstnm']}' does not contain enough data")
#         exit(-1)

#     fig, ax = plt.subplots(nfreqs+1)
#     fig.set_figwidth(8*0.8)
#     fig.set_figheight(8*0.8)
#     ax[0].plot(sac.times, sac.data, label='Unfiltered')
#     ax[0].scatter([tmin, tmax],[0, 0], c='red')
#     ax[0].legend(loc='upper right')
#     ax[0].set_xticks([],[])
#     ax[0].set_yticks([],[])
#     for i in range(nfreqs):
#         filt = gaussian_filters[i]
#         filtered = calculate_narrow_band_if_good(sac, filt, freqs[i], mingroupv, maxgroupv, 1)
#         ax[i+1].plot(sac.times, filtered, label=f'{periods[i]} s')
#         ax[i+1].legend(loc='upper right')
#         ax[i+1].set_xticks([],[])
#         ax[i+1].set_yticks([],[])
#     plt.tight_layout()
#     plt.show()



def test_tdomain_tapering():
    tmin = sac.headers["dist"]/maxgroupv + sac.headers["o"];
    tmax = sac.headers["dist"]/mingroupv + sac.headers["o"];
    taper_width = 0.1 * (tmax - tmin) # 10%
    filt = design_tdomain_taper(sac, tmin, tmax, taper_width)
    filt_data = apply_filter_tdomain(sac.data, filt)
    fig, ax = plt.subplots(3)
    fig.set_figwidth(8*0.8)
    fig.set_figheight(5*0.8)
    ax[0].plot(sac.times, sac.data, label='Original', zorder=1)
    ax[0].scatter([tmin, tmax], [0,0], c='red', zorder=2)
    ax[0].legend(loc='best')
    ax[1].plot(sac.times, filt, label='Hanning taper: %.0f s' %(taper_width))
    ax[1].legend(loc='best')
    ax[2].plot(sac.times, filt_data, label='Filtered Data')
    ax[2].legend(loc='best')
    plt.tight_layout()
    plt.show()


def test_multichannel_isolation_filtering():
    tmin = sac.headers["dist"]/maxgroupv + sac.headers["o"];
    tmax = sac.headers["dist"]/mingroupv + sac.headers["o"];

    taper_width = 0.1 * (tmax - tmin) # 10%
    hanning_filt = design_tdomain_taper(sac, tmin, tmax, taper_width)
    fax, gaussian_filters = design_fdomain_gaussians(freqs,
                            sac.headers['delta'],
                            len(sac.data),
                            filter_width_range)

    fig, ax = plt.subplots(nfreqs+1)
    fig.set_figwidth(6*0.8)
    fig.set_figheight(8*0.8)
    ax[0].plot(sac.times, sac.data, label='Unfiltered')
    ax[0].legend(loc='upper right')
    ax[0].set_xticks([],[])
    ax[0].set_yticks([],[])
    for i in range(nfreqs):
        filtered = apply_filter_fdomain(sac.data, gaussian_filters[i])
        filtered = apply_filter_tdomain(filtered, hanning_filt)
        ax[i+1].plot(sac.times, filtered, label=f'{periods[i]} s')
        ax[i+1].legend(loc='upper right')
        ax[i+1].set_xticks([],[])
        ax[i+1].set_yticks([],[])
    plt.tight_layout()
    plt.show()


def test_sac_class():
    import obspy
    st = obspy.read(test_sacfile)
    tr = st[0]

    sac_data = tr.data
    sac_times = tr.times()
    sac_headers = {}
    for key in tr.stats.sac.keys():
        sac_headers[key] = tr.stats.sac[key]

    sac_obj = Sac()
    sac_obj.append(sac_data, sac_times, sac_headers)
    sac_obj.write('test.sac')

    sac_obj = Sac()
    sac_obj.read(test_sacfile)
    sac_obj.write('test2.sac')




    
    


if __name__=="__main__":
    # test_plot_egfs()
    # test_plot_event()
    # test_narrow_band_gaussian_filter_design()
    # test_narrow_band_filtered_seisomgrams()
    # test_tdomain_tapering()
    # test_multichannel_isolation_filtering()
    test_sac_class()




