#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
#                              Object parameters                           #
############################################################################

#   Name of the variable star
name_star = "?"

#   Coordinates - Format:  ra = hh:mm:ss e.g. 19:44:42.8539591894
#                         dec = dd:am:as e.g. +54:49:42.887193554
ra_obj = "??:??:??"
dec_obj = "+??:??:??"

#   Date of the minimum (UTC)
#   "yyyy:mm:ddThh:mm:ss" e.g., "2020-09-18T01:00:00"
transit_time = "?"

#   Period (Algol: p=2.867315d, RZ Cas: p=1.1952499d, TV Cas: p=1.81259d)
period = '?'

############################################################################
#                Additional options: only edit if necessary                #
############################################################################

############################################################################
#   Light curve parameters
#
#   Filter
filter_ = '?'

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir = 'output/'

#   Light curve file name
file_name = f'{output_dir}/tables/light_curve_{filter_}.csv'

#   Binning in days
binning_factor = 0.0001

############################################################################
#                               Libraries                                  #
############################################################################

from ost_photometry.analyze import plot

from astropy.timeseries import TimeSeries

############################################################################
#                                  Main                                    #
############################################################################

if __name__ == '__main__':
    #   Load light curve
    ts = TimeSeries.read(
        file_name,
        format='ascii.csv',
        time_column='time',
        time_format='jd',
    )

    #   Plot light curve over JD
    plot.light_curve_jd(
        ts,
        filter_,
        f"{filter_}_err",
        output_dir,
        name_obj=name_star,
    )

    #   Plot the light curve folded on the period
    plot.light_curve_fold(
        ts,
        filter_,
        f"{filter_}_err",
        output_dir,
        transit_time,
        period,
        binning_factor=binning_factor,
        name_obj=name_star,
    )
