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
period: float | str = '?'

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

#   Optional input path. If None, uses the default file name below.
#   Can be a single CSV file path or a directory containing multiple CSV files.
input_path: str | None = None

#   Light curve file name (default single-file input)
color_suffix: str = f'_{color}' if color not in (None, '') else ''
file_name: str = f'{output_dir}/tables/light_curve_{name_star.replace(" ", "_")}_{filter_}{color_suffix}.csv'

#   Binning in days
binning_factor: float = 0.0001

############################################################################
#                               Libraries                                  #
############################################################################

from ost_photometry.analyze import plots
import os
from glob import glob
from astropy.table import vstack
import re

from astropy.timeseries import TimeSeries

############################################################################
#                                  Main                                    #
############################################################################

def infer_filter_from_filename(path: str) -> str | None:
    base = os.path.basename(path)
    match = re.search(r'^light_curve_.*?_([A-Za-z]+)(?:_|\.|$)', base)
    return match.group(1) if match else None


def infer_filter_from_timeseries(ts_in: TimeSeries) -> str | None:
    candidate_filters = []
    for column_name in ts_in.colnames:
        if column_name == 'time' or column_name.endswith('_err'):
            continue
        if f"{column_name}_err" in ts_in.colnames:
            candidate_filters.append(column_name)
    if not candidate_filters:
        return None
    if filter_ not in (None, '?') and filter_ in candidate_filters:
        return filter_
    return sorted(candidate_filters)[0]


def resolve_input_files() -> list[str]:
    if input_path is None:
        return [file_name]
    if os.path.isdir(input_path):
        return sorted(glob(os.path.join(input_path, '*.csv')))
    return [input_path]


def plot_single_file(csv_path: str) -> None:
    ts = TimeSeries.read(
        csv_path,
        format='ascii.csv',
        time_column='time',
        time_format='jd',
    )
    ts.sort('time')

    plots.light_curve_jd(
        ts,
        filter_,
        f"{filter_}_err",
        output_dir,
        name_object=name_star,
    )

    if (transit_time is not None and transit_time != '?'
            and period is not None and period != '?' and period > 0.0):
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


def plot_directory(files: list[str]) -> None:
    filter_to_ts_list: dict[str, list[TimeSeries]] = {}
    for _path in files:
        ts_i = TimeSeries.read(
            _path,
            format='ascii.csv',
            time_column='time',
            time_format='jd',
        )

        inferred = infer_filter_from_filename(_path)
        if inferred is None or (inferred not in ts_i.colnames or f"{inferred}_err" not in ts_i.colnames):
            inferred = infer_filter_from_timeseries(ts_i)

        if inferred is None:
            raise ValueError(f"Could not infer filter for file: {_path}")

        filter_to_ts_list.setdefault(inferred, []).append(ts_i)

    for current_filter, ts_parts in sorted(filter_to_ts_list.items()):
        ts_combined = ts_parts[0] if len(ts_parts) == 1 else vstack(ts_parts, metadata_conflicts='silent')
        ts_combined.sort('time')

        plots.light_curve_jd(
            ts_combined,
            current_filter,
            f"{current_filter}_err",
            output_dir,
            name_object=name_star,
        )


        if (transit_time is not None and transit_time != '?'
            and period is not None and period != '?' and period > 0.0):
            plots.light_curve_fold(
                ts_combined,
                current_filter,
                f"{current_filter}_err",
                output_dir,
                transit_time,
                period,
                binning_factor=binning_factor,
                name_object=name_star,
            )


def main() -> None:
    input_files = resolve_input_files()
    if not input_files:
        raise FileNotFoundError('No CSV files found to load the light curve.')

    if len(input_files) > 1 or (input_path is not None and os.path.isdir(input_path)):
        plot_directory(input_files)
    else:
        plot_single_file(input_files[0])


if __name__ == '__main__':
    main()
