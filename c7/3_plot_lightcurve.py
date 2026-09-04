#! /usr/bin/python
# -*- coding: utf-8 -*-

"""Plot light curves from the pipeline table ``output/tables/light_curves.ecsv``.

To overlay several nights, use ``5_compare_nights.py``.
"""

############################################################################
#                              Object parameters                           #
############################################################################

#   Name of the variable star
name_star: str = "?"

#   Coordinates - Format:  ra = hh:mm:ss e.g. 19:44:42.8539591894
#                         dec = dd:am:as e.g. +54:49:42.887193554
#               - used only if object_id is None (nearest source in the table)
ra_obj: str = "??:??:??"
dec_obj: str = "+??:??:??"

object_id: int | None = None

#   Date of the minimum (UTC)
#   "yyyy:mm:ddThh:mm:ss" e.g., "2020-09-18T01:00:00"
transit_time: str = "?"

#   Period (Algol: p=2.867315d, RZ Cas: p=1.1952499d, TV Cas: p=1.81259d)
period: float | str = "?"

############################################################################
#                Additional options: only edit if necessary                #
############################################################################

############################################################################
#   Light curve parameters
#
#   Filter
filter_: str = "?"

#   Colour in the table (e.g. "B-V") if the pipeline wrote colour rows; or "" / None
color: str | None = None

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir: str = "output"
# Default: <output_dir>/tables/light_curves.ecsv

#   Optional input path. If None, auto-detect under ``output_dir/tables/``.
input_path: str | None = None


#   Binning in days
binning_factor: float = 0.0001
binning_factor: float | None = None

#   Output file type
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
    object_id_from_epoch_native_sky,
    plot_from_light_curves_table,
)

############################################################################
#                                  Main                                    #
############################################################################


def _resolve_table_path() -> str:
    if input_path is not None:
        if os.path.isfile(input_path):
            return input_path
        raise FileNotFoundError(f"Light-curve table not found: {input_path!r}")
    path = os.path.join(output_dir, "tables", "light_curves.ecsv")
    if os.path.isfile(path):
        return path
    raise FileNotFoundError(
        f"No {path!r}. Run ``2_obtain_flux.py`` with skip_light_curve=False first."
    )


def _optional_sky() -> SkyCoord | None:
    ra_s, dec_s = ra_obj.strip(), dec_obj.strip()
    if (not ra_s or "?" in ra_s) and (not dec_s or "?" in dec_s):
        return None
    if "?" in ra_s or "?" in dec_s or not ra_s or not dec_s:
        raise ValueError("Provide both sexagesimal RA and Dec, or leave both as placeholders.")
    return SkyCoord(ra_s, dec_s, unit=(u.hourangle, u.deg), frame="icrs")


def _resolve_id(tbl: Table) -> int:
    if object_id is not None:
        return int(object_id)
    sky = _optional_sky()
    if sky is not None:
        return object_id_from_epoch_native_sky(tbl, sky)
    names = np.asarray(tbl["object_name"]).astype(str) if "object_name" in tbl.colnames else None
    if names is not None and name_star not in (None, "?", ""):
        want = name_star.strip()
        for sid, nm in zip(
            np.asarray(tbl["id"]).astype(int),
            names,
            strict=False,
        ):
            if str(nm).strip() == want or str(nm).strip().replace(" ", "_") == want.replace(
                " ", "_"
            ):
                return int(sid)
    raise ValueError(
        "Set object_id, sky coordinates, or name_star matching object_name in light_curves.ecsv."
    )


def _filters_to_plot(tbl: Table) -> list[str]:
    available = sorted(set(np.asarray(tbl["filter"]).astype(str)))
    out: list[str] = []
    if filter_ not in (None, "?", ""):
        if filter_ not in available:
            raise ValueError(
                f"Filter {filter_!r} not in light_curves.ecsv columns {available}"
            )
        out.append(filter_)
    if color not in (None, "", "?-?", "?"):
        if color in available:
            out.append(str(color))
        else:
            print(f"Note: colour {color!r} not in table (have {available}).", file=sys.stderr)
    if not out:
        # Prefer photometric bands over colour strings
        bands = [f for f in available if "-" not in f]
        out = bands or available
    return list(dict.fromkeys(out))


def main() -> None:
    path = _resolve_table_path()
    tbl = Table.read(path)
    sid = _resolve_id(tbl)
    per = period
    if per == "?" or per is None:
        per = None
    else:
        try:
            per = float(per)
        except (TypeError, ValueError):
            per = None
    tt = None if transit_time in (None, "?") else transit_time
    name = name_star if name_star not in (None, "?") else str(sid)

    for filt in _filters_to_plot(tbl):
        plot_from_light_curves_table(
            tbl,
            sid,
            filt,
            output_dir,
            name_object=name,
            file_type=file_type,
            transit_time=tt,
            period=per,
            binning_factor=binning_factor,
        )


if __name__ == "__main__":
    main()
