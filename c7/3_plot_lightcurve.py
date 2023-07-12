#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

####################  Object parameters  ####################

###
#   Name of the variable star
#
namestar     = "?"


###
#   Coordinates - Format:  ra = hh:mm:ss e.g. 19:44:42.8539591894
#                         dec = dd:am:as e.g. +54:49:42.887193554
#
ra_obj       = "??:??:??"
dec_obj      = "+??:??:??"


###
#   Date of the minimum (UTC)
#   "yyyy:mm:ddThh:mm:ss" e.g., "2020-09-18T01:00:00"
#
transit_time = "?"


###
#   Period (Algol: p=2.867315d, RZ Cas: p=1.1952499d, TV Cas: p=1.81259d)
#
period       = ?


############################################################################
####             Additional options: only edit if necessary             ####
############################################################################

##################  Light curve parameters  #################

###
#   Filter
#
filt = '?'


###
#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
#
outdir = 'output/'


###
#   Light curve file name
#
fname = outdir+'/tables/light_curce_'+filt+'.csv'


###
#   Light curve options
#
#   Binning in days
tbin = 0.0001


############################################################################
####                            Libraries                               ####
############################################################################

import sys

from ost_photometry.style import bcolors

from ost_photometry.analyze import plot

from astropy.timeseries import TimeSeries
from astropy.utils.data import get_pkg_data_filename

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Load light curve
    ts = TimeSeries.read(
        fname,
        format='ascii.csv',
        time_column='time',
        time_format='jd',
        )

    #   Plot light curve over JD
    plot.light_curve_jd(ts, filt, filt+'_err', outdir, nameobj=namestar)

    #   Plot the light curve folded on the period
    plot.light_curve_fold(
        ts,
        filt,
        filt+'_err',
        outdir,
        transit_time,
        period,
        binn=tbin,
        nameobj=namestar,
        )
