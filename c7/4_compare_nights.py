#! /usr/bin/python
# -*- coding: utf-8 -*-

"""Compare C7 light curves from different nights.

Reads ``tables/light_curves.ecsv`` from one or more pipeline output
directories (or explicit ECSV paths), matches the variable by id / sky /
name, and writes overlay + per-night panel plots.

A single table that already spans several nights works too: colour is
``night_id`` from JD, not the file name.
"""

############################################################################
#                              Object parameters                           #
############################################################################

#   Name of the variable star
name_star: str = "?"

#   Coordinates — used if object_id is None (nearest source in each table)
#   Format: ra = hh:mm:ss  e.g. 19:44:42.85
#           dec = dd:am:as e.g. +54:49:42.89
ra_obj: str = "??:??:??"
dec_obj: str = "+??:??:??"

object_id: int | None = None

#   Date of the minimum (UTC) — needed for the folded overlay
#   "yyyy-mm-ddThh:mm:ss" e.g. "2020-09-18T01:00:00"
transit_time: str = "?"

#   Period (Algol: 2.867315 d, RZ Cas: 1.1952499 d, TV Cas: 1.81259 d)
period: float | str = "?"

############################################################################
#                Additional options: only edit if necessary                #
############################################################################

#   Filter (must exist in the tables). Leave "?" to use the first band.
filter_: str = "?"

#   Output directory for the comparison PDFs
#   (``<output_dir>/results/lightcurves/``).
output_dir: str = "output_nights"

#   Pipeline output directories (each with tables/light_curves.ecsv) or
#   direct paths to light_curves.ecsv. One entry is enough if that table
#   already contains several nights.
nights: list[str] = [
    "output",
    # "output_night2",
]

#   Subtract / divide the per-night median so shapes line up when the
#   zero point differs between reductions. Default off: keep calibrated mag.
align_nightly_zero: bool = False

file_type: str = "pdf"

############################################################################
#                               Libraries                                  #
############################################################################

import os
import sys

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

from ost_photometry.analyze.post_processing.light_curve import (
    combine_light_curve_night_tables,
    plot_nights_compare_from_table,
    resolve_light_curves_ecsv_path,
)

############################################################################
#                                  Main                                    #
############################################################################


def _optional_sky() -> SkyCoord | None:
    ra_s, dec_s = ra_obj.strip(), dec_obj.strip()
    if (not ra_s or "?" in ra_s) and (not dec_s or "?" in dec_s):
        return None
    if "?" in ra_s or "?" in dec_s or not ra_s or not dec_s:
        raise ValueError(
            "Provide both sexagesimal RA and Dec, or leave both as placeholders."
        )
    return SkyCoord(ra_s, dec_s, unit=(u.hourangle, u.deg), frame="icrs")


def _filters_available(paths: list[str]) -> list[str]:
    found: list[str] = []
    for raw in paths:
        tbl = Table.read(resolve_light_curves_ecsv_path(raw))
        found.extend(np.asarray(tbl["filter"]).astype(str).tolist())
    return sorted(set(found))


def _filters_to_plot(paths: list[str]) -> list[str]:
    available = _filters_available(paths)
    if filter_ not in (None, "?", ""):
        if filter_ not in available:
            raise ValueError(
                f"Filter {filter_!r} not in night tables {available}"
            )
        return [filter_]
    bands = [f for f in available if "-" not in f]
    return bands or available


def main() -> None:
    if not nights:
        raise ValueError("Set ``nights`` to at least one output directory or ECSV.")
    for raw in nights:
        resolve_light_curves_ecsv_path(raw)

    per = period
    if per == "?" or per is None:
        per = None
    else:
        try:
            per = float(per)
        except (TypeError, ValueError):
            per = None
    tt = None if transit_time in (None, "?") else transit_time
    name = name_star if name_star not in (None, "?") else "object"
    sky = _optional_sky()
    name_arg = name_star if name_star not in (None, "?", "") else None

    os.makedirs(output_dir, exist_ok=True)

    for filt in _filters_to_plot(nights):
        tbl = combine_light_curve_night_tables(
            nights,
            filt,
            object_id=object_id,
            coord=sky,
            object_name=name_arg,
        )
        plot_nights_compare_from_table(
            tbl,
            int(tbl["id"][0]),
            filt,
            output_dir,
            name_object=name,
            file_type=file_type,
            transit_time=tt,
            period=per,
            align_nightly_zero=align_nightly_zero,
        )


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as exc:
        print(exc, file=sys.stderr)
        sys.exit(1)
