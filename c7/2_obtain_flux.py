#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
Multi-epoch / multi-filter photometry for variable stars (C7).

Runs ``Observation.run_pipeline`` with ``extraction_mode="multi"`` (directories
per filter), then calibration and ``LightCurveStep`` (CSV under
``<output_dir>/tables/light_curve_*.csv``, plots under ``<output_dir>/results/lightcurves/``).

Calibration via ``PipelineConfig.from_preset`` (``linear_fit_per_image``, ``linear_fit_per_image_extinction``)
or fine-grained ``calibration_strategy`` / ``calibration_grouping`` / ``extinction_mode`` fields.
"""

############################################################################
#             Configuration: modify the file in this section               #
############################################################################

###
#   Name of the variable star
#
name_star: str = "?"

###
#   Coordinates - Format:  ra = hh:mm:ss e.g. 19:44:42.8539591894
#                         dec = dd:am:as e.g. +54:49:42.887193554
#
ra_star: str = "??:??:??"
dec_star: str = "??:??:??"

###
#   Date of the minimum (UTC)
#   "yyyy-mm-ddThh:mm:ss" e.g., "2020-09-18T01:00:00"
#   This parameter is optional. Set it to '?' if unknown.
#
transit_time: str = "?"

###
#   Period (Algol: p=2.867315d, RZ Cas: p=1.1952499d, TV Cas: p=1.81259d)
#   This parameter is optional. Set it to '?' if unknown.
#
period: float | str = '?'


############################################################################
#                Additional options: only edit if necessary                #
#              Add up to 10 variable set for different filter              #
############################################################################

############################################################################
#   Finder options
#
#   Set FWHM -> characterizes the size of the diffraction patterns
#   ``None`` = estimate FWHM per image; set a float to force that value (px).
fwhm: float | None = None

#   Accepted range for automatic FWHM estimation (pixels). Raise the max if
#   seeing is poor and auto-FWHM falls back to the default too often.
fwhm_estimate_min: float = 2.0
fwhm_estimate_max: float = 15.0

#   Detection threshold = this × MAD-RMS of the sky-subtracted frame.
#   Single-epoch C7 frames need ~10 (package default 5 is for stacked N2).
#   Lower this if faint stars vanish; raise it if the finder still tags noise.
multiplier_background_rms: float = 10.0
#   Same for residual re-find in ePSF photometry (only used if method is PSF).
multiplier_background_rms_epsf: float = 10.0

#   Worker processes for multi-image extraction. ``None`` = half the CPUs.
n_cores_multiprocessing: int | None = 8

############################################################################
#   Define filter 1 (e.g., U, B, V, Clear, ...)
#   Clear / C / luminance / white have no catalog ZP. The pipeline then
#   writes a relative flux light curve (continuum ≈ 1) instead of mag_cal_*.
#
filter_1: str = 'V'

############################################################################
#   Path to the images of filter 1
#
path_1: str = './output/V/'

############################################################################
#   Define filter 2 (e.g., U, B,V,...)
#
filter_2: str = 'B'

############################################################################
#   Path to the images of filter 2
#
path_2: str = './output/B/'

############################################################################
#   Define filter 3 (e.g., U, B,V,...)
#
filter_3: str = 'R'

############################################################################
#   Path to the images of filter 3
#
path_3: str = './output/R/'

############################################################################
#   Additional program options
#
#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir: str = 'output/'

#   Aperture or ePSF photometry
photometry_extraction_method: str = 'APER'
# photometry_extraction_method: str = 'PSF'

###
#   Calibration configuration mode: "preset" or "custom"
#
calibration_config_mode: str = "preset"

###
#   Preset (used when calibration_config_mode == "preset")
#
#   ``linear_fit_per_image`` — T/ZP per visit, no extinction (C7 default).
#   ``linear_fit_per_image_extinction`` — same with fitted extinction when airmass varies.
#   ``linear_fit_per_night`` — one T/ZP for the whole night (derive-transform).
#
calibration_preset: str = "linear_fit_per_image"

###
#   Fine-grained calibration (used when calibration_config_mode == "custom")
#
calibration_strategy: str = "linear_fit"      # median_zp | linear_fit
calibration_grouping: str = "per_image"       # per_image | per_night | ensemble | fixed
extinction_mode: str = "none"                 # none | tabulated | from_comparison_stars | from_value_airmass
color_term_fit: str = "auto"                  # always | auto | never
fit_sigma_clip: float = 2.5
derive_transform_from_data: bool = False
zp_subsample_statistic: bool = False
exposure_pairing: str = "jd_nearest"          # jd_nearest | index
exposure_jd_tolerance: float = 0.001
reference_filter: str | None = None

############################################################################
#   Calibration source (possibilities: simbad_vot, UCAC4, GSC2.3, URAT1,
#                                      NOMAD, HMUBV, GSPC2.4, APASS,
#                                      Swift/UVOT, XMM-OM, VRI-NCC)
#
'''
    UCAC4 Catalogue (Zacharias+, 2012):
    B & V data are from the AAVSO Photometric all-sky survey (APASS) DR6
    plus single observation stars kindly provided by A.Henden. For bright
    stars the Bmag and Vmag columns contain the Hipparcos/Tycho Bt and Vt
    mags respectively.
'''
# calibration_source: str = 'UCAC4'


'''
The Full GSC2.3.2 Catalogue:
    Often not Johnson B & V magnitudes, but similar passbands. Only those
    with codes 3 and 4 are Johnson magnitudes.
'''
# calibration_source: str = 'GSC2.3'


'''
    URAT1 Catalog (Zacharias+ 2015)
    B and V are also from the APASS survey.
'''
# calibration_source: str = 'URAT1


'''
    NOMAD Catalog (Zacharias+ 2005):
    V: Photometric magnitude in Optical V band between 500 and 600 nm
    B: Photometric magnitude in Optical B band between 400 and 500 nm
'''
# calibration_source: str = 'NOMAD'


'''
    Homogeneous Means in the UBV System (Mermilliod 1991):
    Johnson V-band & Johnson B-band magnitude (only B-V given)
'''
# calibration_source: str = 'HMUBV'


'''
    Guide Star Photometric Catalog V2.4 (Bucciarelli+ 2001):
    Johnson B,V,R-band magnitude
'''
# calibration_source: str = 'GSPC2.4'


'''
    AAVSO Photo. All Sky Surv. DR9(Henden+,2016):
    Johnson V-band & Johnson B-band magnitude
'''
calibration_source: str = 'APASS'

'''
    Swift/UVOT Serendipitous Source Catalog (Yershov, 2015):
    AB magnitudes: U-AB, B-AB, V-AB
'''
# calibration_source: str = 'Swift/UVOT'


'''
    XMM-OM Serendipitous Source Survey Catalogue (XMM-SUSS5.0)
    (Page+, 2021):
    AB magnitudes: UmAB, BmAB, VmAB
'''
# calibration_source: str = 'XMM-OM'


'''
    Optical-UV-IR survey of North Celestial Cap (Gorbikov+, 2014):
    Johnson VRI magnitudes
'''
# calibration_source: str = 'VRI-NCC'


'''
    The USNO-B1.0 Catalog (Monet+ 2003):
    BRI magnitudes
'''
# calibration_source: str = 'USNO-B1.0'


#   Magnitude limit of the calibration stars
magnitude_range: tuple[float, float] = (12., 15.)

#   On-sky match to the calibration catalog (arcsec). Independent of
#   inter-filter ``separation_limit`` (default 1"). Tighten in crowded fields
#   to cut the catalog-crossmatch tail.
calibration_match_radius_arcsec: float = 1.0

############################################################################
#   Aperture options
#
#   If True, APER radii are factor × image FWHM (pixels). The absolute
#   radii below are then unused for extraction. Limiting-magnitude
#   ``aperture_radius`` still uses ``radius_aperture``. Typical factor
#   1.5–2.5; annulus factors must stay larger than the aperture factor.
aperture_scale_with_fwhm: bool = False
aperture_fwhm_factor: float = 2.0
inner_annulus_fwhm_factor: float = 2.8
outer_annulus_fwhm_factor: float = 4.0

#   Extraction radius stars in arcsec or pixel (used when the flag is False)
radius_aperture: float = 5.

#   Extraction radius background (inner and outer radii) in arcsec or pixel
inner_annulus_radius: float = 7.
outer_annulus_radius: float = 10.

#   Unit
# r_unit: str = 'pixel'
radii_unit: str = 'arcsec'

############################################################################
#   Correlation options
#
#   ID of the reference image (``0`` or ``"auto"`` = most detections / sharpest)
reference_image_index: int | str = "auto"

#   Maximal sky separation between two objects in arcsec (inter-filter and
#   unaligned intra-filter). Tightening this does not hide HAT-P-20: OOI
#   identification uses at least 5" unless you set ooi_separation_limit.
separation_limit: float = 5.

#   Sparse tracks: keep stars that miss some frames (long C7 series).
require_complete_intersection: bool = False
#   Keep a track if it is detected on at least this fraction of frames.
min_detection_fraction: float = 0.3
#   Miss-count floor; combined with the fraction so a small number does not
#   force near-completeness on a long series.
n_allowed_non_detections_object: int = 5

#   Independent ASTAP solutions on aligned frames make the same pixel look
#   like a different sky position and mix track IDs. Keep False when
#   ``1_reduce_images.py`` used shift_all=True (``aa_true`` or ``wcs``).
#   True only for native (unwarped) pixels.
wcs_solve_all_images: bool = False

#   Intra-filter tracks: ``auto`` uses pixel matching when the series is
#   registered (``aa_true`` or ``wcs`` with shift_all). ``sky`` is for
#   unaligned native extract.
correlation_coordinates: str = "auto"
# correlation_coordinates: str = "pixel"
# correlation_coordinates: str = "sky"

#   ``to_reference`` with pixel matching keeps identity on the reference
#   grid. ``sequential`` is a fallback when frames are not registered.
correlation_link_mode: str = "to_reference"
# correlation_link_mode: str = "sequential"

############################################################################
#   Light curve options
#
#   Phase-bin width for folded curves (0–1), or None for no binning
# binning_factor: float | None = 0.0001
# binning_factor: float | None = 0.0002
binning_factor: float | None = None

############################################################################
#                               Libraries                                  #
############################################################################

import time

import warnings

warnings.filterwarnings('ignore')

from ost_photometry import style
from ost_photometry.analyze import analyze
from ost_photometry.analyze.pipeline import PipelineConfig

import astropy.units as u

############################################################################
#                                  Main                                    #
############################################################################

if __name__ == '__main__':
    #   Set start time
    start_time = time.time()

    #   Prepare variable lists and dictionaries from the individual
    #   definitions above
    filter_list: list[str] = []
    image_paths: dict[str, str] = {}
    fwhm_object_psf: dict[str, float] = {}
    for i in range(0, 10):
        if 'filter_' + str(i) in locals():
            filter_list.append(locals()['filter_' + str(i)])
            image_paths[locals()['filter_' + str(i)]] = locals()['path_' + str(i)]
            fwhm_object_psf[locals()['filter_' + str(i)]] = fwhm

    fwhm_for_pipeline: dict[str, float] | None = None
    if any(v is not None for v in fwhm_object_psf.values()):
        fwhm_for_pipeline = {
            k: float(v) for k, v in fwhm_object_psf.items() if v is not None
        }

    ###
    #   Initialize observation and run full pipeline (incl. light curves)
    #
    observation = analyze.Observation(
        ra_objects=[ra_star],
        dec_objects=[dec_star],
        object_names=[name_star],
        ra_unit=u.hourangle,
        dec_unit=u.deg,
        transit_times=[transit_time],
        periods=[period],
    )

    _color_idx = None
    if len(filter_list) >= 2:
        f_a, f_b = filter_list[0], filter_list[1]
        _color_idx = {f_b: (f_a, f_b)}

    _mode = calibration_config_mode.strip().lower()
    if _mode not in ("preset", "custom"):
        raise ValueError(
            "calibration_config_mode must be 'preset' or 'custom', "
            f"got {calibration_config_mode!r}"
        )

    _shared_pipeline_kw = dict(
        fwhm_object_psf=fwhm_for_pipeline,
        fwhm_estimate_min=fwhm_estimate_min,
        fwhm_estimate_max=fwhm_estimate_max,
        n_cores_multiprocessing=n_cores_multiprocessing,
        multiplier_background_rms=multiplier_background_rms,
        multiplier_background_rms_epsf=multiplier_background_rms_epsf,
        photometry_extraction_method=photometry_extraction_method,
        aperture_scale_with_fwhm=aperture_scale_with_fwhm,
        aperture_fwhm_factor=aperture_fwhm_factor,
        inner_annulus_fwhm_factor=inner_annulus_fwhm_factor,
        outer_annulus_fwhm_factor=outer_annulus_fwhm_factor,
        radius_aperture=radius_aperture,
        inner_annulus_radius=inner_annulus_radius,
        outer_annulus_radius=outer_annulus_radius,
        radii_unit=radii_unit,
        reference_image_index=reference_image_index,
        n_allowed_non_detections_object=n_allowed_non_detections_object,
        require_complete_intersection=require_complete_intersection,
        min_detection_fraction=min_detection_fraction,
        wcs_solve_all_images=wcs_solve_all_images,
        correlation_link_mode=correlation_link_mode,
        correlation_coordinates=correlation_coordinates,
        separation_limit=separation_limit * u.arcsec,
        calibration_source=calibration_source,
        calibration_catalog_mag_range=magnitude_range,
        calibration_match_radius=calibration_match_radius_arcsec * u.arcsec,
        color_indices=_color_idx,
        light_curve_color=(
            f"{filter_list[0]}-{filter_list[1]}" if len(filter_list) >= 2 else None
        ),
        aperture_radius=radius_aperture,
        # Magnitude output stays in the calibrated catalog system
        # (APASS → Johnson/Vega). To convert: convert_magnitudes=True plus
        # e.g. output_filter_set="sdss" (SDSS+Vega is rejected).
        skip_light_curve=False,
        light_curve_binning_factor=binning_factor,
        plot_light_curve_calibration_objects=True,
        skip_derive_limiting_magnitude=True,
        # Cap individual inter-filter pair PDFs (None=all, 0=overview only):
        # diagnostic_plots__correlation_inter_filter_max_pair_plots=25,
        # Cap per-epoch calibration QC PDFs (None=all, 0=skip; default 25):
        # diagnostic_plots__calibration_max_epoch_plots=25,
    )

    if _mode == "preset":
        config = PipelineConfig.from_preset(
            calibration_preset.strip(),
            overrides=_shared_pipeline_kw,
        )
    else:
        _strategies = frozenset({"median_zp", "linear_fit"})
        _groupings = frozenset({"per_image", "per_night", "ensemble", "fixed"})
        _extinction = frozenset({
            "none",
            "tabulated",
            "from_comparison_stars",
            "from_value_airmass",
        })
        _color_term_fits = frozenset({"always", "auto", "never"})
        _pairing = frozenset({"jd_nearest", "index"})
        _cs = calibration_strategy.strip().lower()
        _cg = calibration_grouping.strip().lower()
        _em = extinction_mode.strip().lower()
        _ctf = color_term_fit.strip().lower()
        _ep = exposure_pairing.strip().lower()
        if _cs not in _strategies:
            raise ValueError(f"calibration_strategy must be one of {sorted(_strategies)}")
        if _cg not in _groupings:
            raise ValueError(f"calibration_grouping must be one of {sorted(_groupings)}")
        if _em not in _extinction:
            raise ValueError(f"extinction_mode must be one of {sorted(_extinction)}")
        if _ctf not in _color_term_fits:
            raise ValueError(
                f"color_term_fit must be one of {sorted(_color_term_fits)}"
            )
        if _ep not in _pairing:
            raise ValueError(f"exposure_pairing must be one of {sorted(_pairing)}")
        config = PipelineConfig(
            calibration_strategy=_cs,
            calibration_grouping=_cg,
            extinction_mode=_em,
            color_term_fit=_ctf,
            fit_sigma_clip=fit_sigma_clip,
            derive_transform_from_data=derive_transform_from_data,
            zp_subsample_statistic=zp_subsample_statistic,
            exposure_pairing=_ep,
            exposure_jd_tolerance=exposure_jd_tolerance,
            reference_filter=reference_filter,
            **_shared_pipeline_kw,
        )

    observation.run_pipeline(
        filter_list,
        image_paths=image_paths,
        output_dir=output_dir,
        config=config,
        extraction_mode="multi",
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.))
