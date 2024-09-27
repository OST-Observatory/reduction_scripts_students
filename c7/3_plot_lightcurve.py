#! /usr/bin/python
# -*- coding: utf-8 -*-

############################################################################
#                              Object parameters                           #
############################################################################

#   Name of the variable star
name_star: str = "?"

#   Coordinates - Format:  ra = hh:mm:ss e.g. 19:44:42.8539591894
#                         dec = dd:am:as e.g. +54:49:42.887193554
ra_obj: str = "??:??:??"
dec_obj: str = "+??:??:??"

#   Date of the minimum (UTC)
#   "yyyy:mm:ddThh:mm:ss" e.g., "2020-09-18T01:00:00"
transit_time: str = "?"

#   Period (Algol: p=2.867315d, RZ Cas: p=1.1952499d, TV Cas: p=1.81259d)
period: float = '?'

############################################################################
#                Additional options: only edit if necessary                #
############################################################################

############################################################################
#   Light curve parameters
#
#   Filter
filter_: str = '?'

#   Color - Format "filter_1-filter_2" such as "B-V", can also be '' or None
color: str = '?-?'

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir: str = 'output'

#   Light curve file name
if color is not None or color != '':
    color = f'_{color}'
file_name: str = f'{output_dir}/tables/light_curve_{name_star.replace(" ", "_")}_{filter_}{color}.csv'

#   Binning in days
binning_factor: float = 0.0001

############################################################################
#                               Libraries                                  #
############################################################################

from ost_photometry.analyze import plots

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
    plots.light_curve_jd(
        ts,
        filter_,
        f"{filter_}_err",
        output_dir,
        name_object=name_star,
    )

    #   Plot the light curve folded on the period
    plots.light_curve_fold(
        ts,
        filter_,
        f"{filter_}_err",
        output_dir,
        transit_time,
        period,
        binning_factor=binning_factor,
        name_object=name_star,
    )
