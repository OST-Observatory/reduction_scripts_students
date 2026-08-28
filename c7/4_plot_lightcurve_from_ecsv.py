#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
Light curves from epoch-native ``calibrated_magnitudes_<APER|PSF>_<F1>-<F2>.ecsv``.

Uses :func:`ost_photometry.analyze.post_processing.light_curve.plot_light_curve_from_epoch_native_ecsv`
(same plotting path as the pipeline ``LightCurveStep``).

JD requirements
---------------
Current pipeline writes an ``observation_jd`` column when saving calibrated ECSV.
Older files need ``--epoch-meta-json`` (same structure as
``context.calibration_epoch_meta``); export once e.g. with
``save_calibration_epoch_meta_json`` from ``ost_photometry.analyze.post_processing``.

Source selection
----------------
Use ``object_id`` (or ``--object-id``), or sky position for nearest ``id``:

- ICRS **degrees**: ``object_ra_deg`` / ``object_dec_deg`` or ``--ra-deg`` /
  ``--dec-deg``
- **Sexagesimal** (same style as ``3_plot_lightcurve.py``): ``object_ra`` /
  ``object_dec`` as ``hh:mm:ss`` and ``dd:am:as``, or ``--ra-hms`` /
  ``--dec-dms``

Limit: ``match_max_sep_arcsec`` / ``--max-match-sep-arcsec``. If both degree
and sexagesimal pairs are set in the config, **degrees take precedence**.

Pipeline-style flux normalization
---------------------------------
Set ``pipeline_flux_normalization`` / ``--pipeline-flux-normalization`` with
``quantity="flux"`` to apply the same quasi flux calibration and per-object
normalization as the pipeline ``LightCurveStep`` flux fallback (without loading
``ImageSeries`` — the ECSV is folded into an epoch×source flux matrix).
"""

############################################################################
#                               Libraries                                  #
############################################################################

import argparse
import os
import sys
from typing import Any, Dict, Optional

import astropy.units as u
from astropy.coordinates import SkyCoord

from ost_photometry.analyze.post_processing.light_curve import (
    plot_light_curve_from_epoch_native_ecsv,
)

############################################################################
#                         Configuration (edit here)                        #
############################################################################

# Epoch-native table (under output_dir/tables/ if relative path)
ecsv_path: str = "output/tables/calibrated_magnitudes_APER_B-V.ecsv"

output_dir: str = "output"

# Correlated source id (ECSV ``id`` column), or use sky position below instead
object_id: int = 0

# Optional ICRS degrees: if both set, nearest table source is used (median ra/dec per id)
object_ra_deg: Optional[float] = None
object_dec_deg: Optional[float] = None

# Optional sexagesimal (overridden by object_ra_deg/object_dec_deg if both floats set)
# Format: RA = hh:mm:ss e.g. 19:44:42.85 ; Dec = dd:am:as e.g. +54:49:42.88
object_ra: str = "??:??:??"
object_dec: str = "+??:??:??"

match_max_sep_arcsec: float = 5.0

name_star: str = "my_star"

# Band for the light curve (must match a column prefix, e.g. mag_cal_V, flux_inst_V)
filter_: str = "V"

# "magnitude" -> mag_cal_* / mag_inst_* ; "flux" -> flux_inst_* / flux_*
quantity: str = "magnitude"

# "auto" | "transformed" | "simple" — only affects magnitude mode when both
# epoch_k and epoch_k_simple rows exist
calibration_rows: str = "auto"

# Optional: JSON from save_calibration_epoch_meta_json (if ECSV has no observation_jd)
epoch_meta_json: str | None = None

transit_time: str = "?"
period: float | str = "?"
binning_factor: float | None = None

# Like LightCurveStep ``_run_flux_fallback_for_filter`` (needs quantity="flux")
pipeline_flux_normalization: bool = False
distribution_samples: int = 1000

############################################################################
#                               Functions                                  #
############################################################################

def _optional_icrs_deg_from_sexagesimal(
    ra_s: str, dec_s: str
) -> Optional[tuple[float, float]]:
    """
    Parse RA (hours:min:sec) and Dec (deg:am:as) to ICRS degrees.

    Returns None if both strings look unset (placeholders). Raises ValueError
    if only one of the pair is set.
    """
    ra_s, dec_s = ra_s.strip(), dec_s.strip()
    ra_ph = (not ra_s) or ("?" in ra_s)
    dec_ph = (not dec_s) or ("?" in dec_s)
    if ra_ph and dec_ph:
        return None
    if ra_ph or dec_ph:
        raise ValueError(
            "Provide both sexagesimal RA and Dec (or leave placeholders in both)."
        )
    c = SkyCoord(ra_s, dec_s, unit=(u.hourangle, u.deg), frame="icrs")
    return float(c.ra.deg), float(c.dec.deg)


def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot light curve from calibrated_magnitudes_*.ecsv"
    )
    p.add_argument(
        "--ecsv",
        default=None,
        help="Path to calibrated epoch-native ECSV (default: config ecsv_path)",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        help="Output root (plots under results/lightcurves/, CSV under tables/)",
    )
    p.add_argument("--filter", default=None, help="Filter name, e.g. V")
    p.add_argument("--object-id", type=int, default=None)
    p.add_argument(
        "--ra-deg",
        type=float,
        default=None,
        help="ICRS RA in degrees; use with --dec-deg to pick source by sky position",
    )
    p.add_argument(
        "--dec-deg",
        type=float,
        default=None,
        help="ICRS Dec in degrees; use with --ra-deg",
    )
    p.add_argument(
        "--ra-hms",
        default=None,
        metavar="H:M:S",
        help='RA as hh:mm:ss (hours); use with --dec-dms, e.g. "19:44:42.85"',
    )
    p.add_argument(
        "--dec-dms",
        default=None,
        metavar="D:M:S",
        help='Dec as dd:am:as; use with --ra-hms, e.g. "+54:49:42.88"',
    )
    p.add_argument(
        "--max-match-sep-arcsec",
        type=float,
        default=None,
        help="Max on-sky separation for RA/Dec match (default: config or 5)",
    )
    p.add_argument("--object-name", default=None)
    p.add_argument(
        "--quantity",
        choices=("magnitude", "flux"),
        default=None,
        help="Y-axis: calibrated magnitude or extracted flux",
    )
    p.add_argument(
        "--calibration",
        choices=("auto", "transformed", "simple"),
        default=None,
        dest="calibration_rows",
        help="Which epoch_id rows for magnitudes (when simple+transformed both exist)",
    )
    p.add_argument(
        "--epoch-meta-json",
        default=None,
        help="Optional JSON for JD lookup if ECSV has no observation_jd column",
    )
    p.add_argument("--binning-factor", type=float, default=None)
    p.add_argument("--transit-time", default=None)
    p.add_argument("--period", default=None)
    p.add_argument(
        "--pipeline-flux-normalization",
        action="store_true",
        help="Quasi flux calibration + per-object median norm (pipeline flux fallback); requires --quantity flux",
    )
    p.add_argument(
        "--distribution-samples",
        type=int,
        default=None,
        help="Samples for flux uncertainty distributions (default: config or 1000)",
    )
    return p

############################################################################
#                                  Main                                    #
############################################################################

def main(argv: list[str] | None = None) -> None:
    args = _build_arg_parser().parse_args(argv)

    path = args.ecsv or ecsv_path
    out = args.output_dir or output_dir
    filt = args.filter or filter_
    ra_cli = args.ra_deg
    dec_cli = args.dec_deg
    ra_hms_cli = args.ra_hms
    dec_dms_cli = args.dec_dms
    ra_cfg = object_ra_deg
    dec_cfg = object_dec_deg
    if (ra_cli is not None) ^ (dec_cli is not None):
        print(
            "Error: provide both --ra-deg and --dec-deg for sky matching.",
            file=sys.stderr,
        )
        sys.exit(1)
    if (ra_hms_cli is not None) ^ (dec_dms_cli is not None):
        print(
            "Error: provide both --ra-hms and --dec-dms for sexagesimal sky matching.",
            file=sys.stderr,
        )
        sys.exit(1)
    if (ra_cfg is not None) ^ (dec_cfg is not None):
        print(
            "Error: set both object_ra_deg and object_dec_deg in config, or neither.",
            file=sys.stderr,
        )
        sys.exit(1)
    max_sep = (
        args.max_match_sep_arcsec
        if args.max_match_sep_arcsec is not None
        else match_max_sep_arcsec
    )
    oid = args.object_id if args.object_id is not None else object_id
    oname = args.object_name or name_star
    qty = args.quantity or quantity
    calrows = args.calibration_rows or calibration_rows
    meta_json = args.epoch_meta_json or epoch_meta_json
    binning = (
        args.binning_factor if args.binning_factor is not None else binning_factor
    )
    tt = args.transit_time if args.transit_time is not None else transit_time
    per = args.period if args.period is not None else period
    pipe_norm = bool(args.pipeline_flux_normalization or pipeline_flux_normalization)
    dist_samples = (
        args.distribution_samples
        if args.distribution_samples is not None
        else distribution_samples
    )

    if pipe_norm and qty != "flux":
        print(
            "Error: pipeline_flux_normalization requires quantity='flux' "
            "(set quantity in config or use --quantity flux).",
            file=sys.stderr,
        )
        sys.exit(1)

    if not os.path.isfile(path):
        # Try relative to output_dir
        alt = os.path.join(out, "tables", os.path.basename(path))
        if os.path.isfile(alt):
            path = alt
        else:
            print(f"Error: ECSV not found: {path!r}", file=sys.stderr)
            sys.exit(1)

    try:
        sex_cfg = _optional_icrs_deg_from_sexagesimal(object_ra, object_dec)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(
            f"Error: could not parse object_ra/object_dec in config ({e!s}).",
            file=sys.stderr,
        )
        sys.exit(1)

    sex_cli_deg: Optional[tuple[float, float]] = None
    if ra_hms_cli is not None and dec_dms_cli is not None:
        try:
            sex_cli_deg = _optional_icrs_deg_from_sexagesimal(
                ra_hms_cli, dec_dms_cli
            )
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(
                f"Error: could not parse --ra-hms/--dec-dms ({e!s}).",
                file=sys.stderr,
            )
            sys.exit(1)
        if sex_cli_deg is None:
            print(
                "Error: --ra-hms and --dec-dms must be valid coordinates.",
                file=sys.stderr,
            )
            sys.exit(1)

    kw: Dict[str, Any] = dict(
        filter_=filt,
        object_name=oname.replace(" ", "_"),
        quantity=qty,
        calibration_rows=calrows,
        epoch_meta_json=meta_json,
        binning_factor=binning,
        transit_time=tt,
        period=per,
        max_match_sep_arcsec=max_sep,
        pipeline_flux_normalization=pipe_norm,
        distribution_samples=int(dist_samples),
    )
    if ra_cli is not None and dec_cli is not None:
        kw["object_ra_deg"] = ra_cli
        kw["object_dec_deg"] = dec_cli
        kw["object_id"] = None
    elif sex_cli_deg is not None:
        kw["object_ra_deg"] = sex_cli_deg[0]
        kw["object_dec_deg"] = sex_cli_deg[1]
        kw["object_id"] = None
    elif args.object_id is not None:
        kw["object_id"] = int(args.object_id)
    elif ra_cfg is not None and dec_cfg is not None:
        kw["object_ra_deg"] = ra_cfg
        kw["object_dec_deg"] = dec_cfg
        kw["object_id"] = None
    elif sex_cfg is not None:
        kw["object_ra_deg"] = sex_cfg[0]
        kw["object_dec_deg"] = sex_cfg[1]
        kw["object_id"] = None
    else:
        kw["object_id"] = int(oid)

    plot_light_curve_from_epoch_native_ecsv(path, out, **kw)


if __name__ == "__main__":
    main()
