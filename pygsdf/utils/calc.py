#!/usr/bin/env python3

# Functions for basic/general calculations

# Coded by: omid.bagherpur@gmail.com
# Update: 16 Oct 2023


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
