#!/usr/bin/env python

import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


usage = f'''\
usage: python {sys.argv[0]} <processed_events_dirs> <output_dir>
'''

#=====Parameters=====#
periods = [20, 25, 32, 40, 50, 60, 80, 100]
comp = 'LHZ'

#===================#

events_dir = []


def freqs_from_periods(periods):
    freqs = []
    for prd in periods:
        freqs.append(1/prd)
    return sorted(freqs)


def build_gaussian_filter(freqs,dt,npts,minwidth,maxwidth):
    '''
    list(float):  freqs
    int: npts
    float: dt, minwidth, maxwidth
    -----
    freqs: central frequencies
    dt: timeseries delta
    npts: number of points
    '''
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


freqs = freqs_from_periods(periods)
fax, gfilt = build_gaussian_filter(freqs, 1, 3000, 0.06, 0.10)

fig = plt.figure(figsize=(10,5))

ax = fig.add_subplot(111)

for i in range(8):
    plt.plot(fax, gfilt[i])

plt.xlim([0, 0.08])

print(freqs)

plt.show()



