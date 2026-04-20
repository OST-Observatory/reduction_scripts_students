#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
Multi-epoch / multi-filter photometry for variable stars (C7).

Runs ``Observation.run_pipeline`` with ``extraction_mode="multi"`` (directories
per filter), then calibration and ``LightCurveStep`` (CSV under
``<output_dir>/tables/light_curve_*.csv``, plots under ``<output_dir>/lightcurve/``).

Calibration branch is selected with ``calibration_module`` in the configuration
section: ``"legacy"`` (classic zero point + transform) or ``"differential"``
(epoch-wise ``PhotometryCalibrator`` / T-ZP).
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
fwhm: float | None = None

############################################################################
#   Define filter 1 (e.g., U, B,V,...)
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
#   Calibration branch of the analysis pipeline
#
#   ``"legacy"`` — zero point (+ optional magnitude transform) via the classic
#   calibration path (``derive_transformation_coefficients`` /
#   ``calculate_zero_point_statistic`` apply).
#
#   ``"differential"`` — differential photometry (``PhotometryCalibrator``): for
#   two filters, a default color index ``{second: (first, second)}`` is set from
#   ``filter_1`` / ``filter_2`` order below; extend ``PipelineConfig`` in code if
#   you need a different color definition.
#
calibration_module: str = 'legacy'

###
#   Differential calibration only (ignored when ``calibration_module == "legacy"``)
#
#   ``differential_coefficient_mode`` — how T / zero point are grouped:
#   ``per_image``, ``per_night``, ``fixed``, ``ensemble`` (see ``PipelineConfig``).
#
#   ``differential_extinction_order`` — extinction correction in the differential
#   path: ``none``, ``first``, ``second``.
#
differential_coefficient_mode: str = 'per_night'
differential_extinction_order: str = 'first'

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

############################################################################
#   Aperture options
#
#   Extraction radius stars in arcsec or pixel
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
#   ID of the reference image
reference_image_index: int = 0

#   Maximal separation between two objects in arcsec
separation_limit: float = 5.

#   Limit for the number of images on which an object is not found.
#   When this limit is reached, the corresponding object is discarded.
n_allowed_non_detections_object: int = 5

############################################################################
#   Light curve options
#
#   Binning in days (set to None to deactivate)
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

    _cal_mod = calibration_module.strip().lower()
    if _cal_mod not in ('legacy', 'differential'):
        raise ValueError(
            "calibration_module must be 'legacy' or 'differential', "
            f'got {calibration_module!r}'
        )

    # Legacy: fit/use classic transformation & ZP stats. Differential: handled
    # inside PhotometryCalibrator; these legacy flags are turned off.
    _legacy_cal = _cal_mod == 'legacy'
    _color_idx = None
    if _cal_mod == 'differential' and len(filter_list) >= 2:
        f_a, f_b = filter_list[0], filter_list[1]
        _color_idx = {f_b: (f_a, f_b)}

    _diff_pipeline_opts: dict[str, str] = {}
    if _cal_mod == 'differential':
        _coeff_modes = frozenset(
            {'per_image', 'per_night', 'fixed', 'ensemble'},
        )
        _ext_orders = frozenset({'none', 'first', 'second'})
        _dcm = differential_coefficient_mode.strip().lower()
        _deo = differential_extinction_order.strip().lower()
        if _dcm not in _coeff_modes:
            raise ValueError(
                'differential_coefficient_mode must be one of '
                f'{sorted(_coeff_modes)}, got {differential_coefficient_mode!r}'
            )
        if _deo not in _ext_orders:
            raise ValueError(
                'differential_extinction_order must be one of '
                f'{sorted(_ext_orders)}, got {differential_extinction_order!r}'
            )
        _diff_pipeline_opts = {
            'differential_coefficient_mode': _dcm,
            'differential_extinction_order': _deo,
        }

    config = PipelineConfig(
        fwhm_object_psf=fwhm_for_pipeline,
        photometry_extraction_method=photometry_extraction_method,
        radius_aperture=radius_aperture,
        inner_annulus_radius=inner_annulus_radius,
        outer_annulus_radius=outer_annulus_radius,
        radii_unit=radii_unit,
        reference_image_index=reference_image_index,
        wcs_method="astrometry",
        max_pixel_between_objects=3,
        ooi_correlation_strategy=1,
        cross_identification_limit=1,
        n_allowed_non_detections_object=n_allowed_non_detections_object,
        separation_limit=separation_limit * u.arcsec,
        calibration_module=_cal_mod,
        calibration_source=calibration_source,
        calibration_catalog_mag_range=magnitude_range,
        apply_transformation=True,
        derive_transformation_coefficients=_legacy_cal,
        calculate_zero_point_statistic=_legacy_cal,
        differential_color_indices=_color_idx,
        **_diff_pipeline_opts,
        aperture_radius=radius_aperture,
        extract_only_circular_region=False,
        identify_cluster_gaia_data=False,
        clean_objs_using_pm=False,
        convert_magnitudes=False,
        skip_light_curve=False,
        light_curve_binning_factor=binning_factor,
        plot_light_curve_objects_of_interest=True,
        plot_light_curve_calibration_objects=True,
        plot_light_curve_all_objects=False,
        skip_extinction_fit=True,
        skip_derive_limiting_magnitude=True,
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
