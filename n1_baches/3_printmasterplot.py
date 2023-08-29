#! /usr/bin/python3
# -*- coding: utf-8 -*-

############################################################################
#             Configuration: modify the file in this section               #
############################################################################

#   Name of file with individual orders
file_with_orders = "master_spectrum_wr.fit"

#   Name of file with merged spectrum
file_with_merged_spectrum = "master_spectrum_wrm.fit"

#   Name of the object
object_name = "star"

###
#   Radial velocity [km/s]
#       The specification of the radial velocity is necessary for the
#       line identification to work correctly (see below).
radial_velocity = 0.

###
#   Line identifications
#
#   Ions for which line markers are to be drawn.
#   Example: ["HI", "FeI", ...]
ions = []

#   Add lines that ar not in the line file
#   Format: {"Element descriptor": [[wavelength, alignment parameter]]}
#           alignment parameter possibilities: "center", "left", "right"
manual_lines = {"Example Element": [[0., "center"]]}

#   Percent the line flux must be lower than the continuum
percentage_line_flux_must_be_below_continuum = 3.

############################################################################
#                Additional options: only edit if necessary                #
############################################################################

###
#   Extract individual orders and merge by means of PYTHON
#   Possibilities: True or False
#
individual_orders = False

###
#   Panel plot version
#   Possibilities: `old` or `default`
#
panel_version = 'default'
# panel_version = 'old'

#   Wavelength range (in Angstrom) for the panels in the plots
#   (panel_version=old)
panel_wave_range = 500

#   Number of panels on each page/plot
#   (panel_version=default)
n_panels_per_page = 5

#   Number of panels into which the spectrum should be split
#   (panel_version=default)
n_panels = 21

###
#   Normalization ?
#   Possibilities: True or False
#
# normalize = False
normalize = True

#   Normalization version
#   Possibilities: `median_max_window` or `specutils_continuum`
#   Default is `median_max_window`
# norm_version = 'specutils_continuum'
norm_version = 'median_max_window'

#   Normalize before merging of the orders
#   Possibilities: True or False
#   If `False` the orders will be merged first and the merged spectrum will
#   be normalized afterward.
#   (norm_version = specutils_continuum)
norm_individual = False
# norm_individual = True

#   Order of the polynomial used to normalize the spectra
#   (norm_version = specutils_continuum)
porder = 5

#   Wavelength window used to estimate continuum points
#   (norm_version = specutils_continuum)
median_window = 61

#   Number of wavelength bins that should be removed from the beginning
#   of the orders to improve flux normalization between the orders
trim_value_start_order = 0
# trim_value_start_order = 400


###
#   Line identifications
#
#   File containing line identifications
# line_file = ""
# line_file = "absorption_lines.dat"
line_file = "/home/pollux/reduction_scripts_students/n1_baches/atomic_lines.tsv"

###
#   Apply barycentric correction?
#
correct_bary = False

###
#   Debug options
#
#   Debug plot: Creates a plot that can be used check order merging
#   Possibilities: True or False
debug_plot = False
# debug_plot = True


############################################################################
#                           Version history                                #
############################################################################

# v2.1 19.01.2023
#  - new panel version (Fabian)
#  - improved normalization (Fabian)
#  - improved spectral line markers (Fabian)
#  - Bug fix in reading individual orders

# v2.0 08.11.2022
#  - complete rewrite

# v1.3 08.11.2020
#  - added option to set maximum plot height for the subplots

# v1.2 30.07.2019
#  - fixed scaling if the flux is negative
#  - improved scaling of the multiplets
#  - improved overall scaling
#  - continuum determination: improved treatment of the boundary values

# v1.1 18.12.2018
#  - Bugfix with ident plot line numbers
#  - switched from PdfPages to pdfunite
#  - individual orders are no also individual files
#  - enhanced output
#  - shifted normalize option out of user input

# v1.0 python3


############################################################################
#                               Libraries                                  #
############################################################################

import os
import copy

import tempfile

from scipy import interpolate, optimize

import math

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.modeling import models, fitting
from astropy.stats import sigma_clip, mad_std
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
import astropy.units as u

from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.ndimage import median_filter, maximum_filter
from scipy.constants import c

from matplotlib.backends.backend_pdf import PdfPages

from specutils import Spectrum1D
from specutils.fitting import fit_generic_continuum
from specutils.spectra import SpectralRegion

from PyPDF2 import PdfMerger


############################################################################
#                               Functions                                  #
############################################################################

class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


plt.rcParams['pdf.fonttype'] = 42

font = {'family': 'serif',
        'color': 'red',
        'weight': 'normal',
        'size': 12,
        }


def consecutive(data, step_size=1):
    """
        Find consecutive elements in a numpy array

        Idee       : stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-in-a-numpy-array
        ----

        Parameters
        ----------
        data       : `numpy.ndarray`
            Data array

        step_size  : `integer`, optional
            Step size between the consecutive elements
    """
    return np.split(data, np.where(np.diff(data) != step_size)[0] + 1)


def bary_correction(fits_file):
    """
        Calculate barycentric velocity correction

        Parameters
        ----------
        fits_file           : `string`
            Name of the file to read

        Returns
        -------
        bvc                 : `float`
            Barycentric velocity correction
    """
    #   Get Header
    header = fits.open(fits_file)[0].header

    #   Extract header information
    try:
        date = header['JD']
        observation_time = Time(date, format='jd', scale='utc')
    except RuntimeError:
        date = header['DATE-OBS']
        observation_time = Time(date, format='fits', scale='utc')

    ra = header['OBJCTRA']
    dec = header['OBJCTDEC']
    observatory_latitude = header.get('SITELAT', 52.4)
    observatory_longitude = header.get('SITELONG', 12.96)
    observatory_height_above_msl = header.get('SITEHIGH', 40.)

    #   Define Location
    observation_site = EarthLocation(
        lat=observatory_latitude,
        lon=observatory_longitude,
        height=observatory_height_above_msl,
    )

    #   Calculate barycentric velocity correction with the help of a
    #   SkyCoord object
    observation_coordinates = SkyCoord(
        ra=ra,
        dec=dec,
        unit=(u.hourangle, u.deg),
        frame='icrs',
        obstime=observation_time,
        location=observation_site
    )

    bvc_value = observation_coordinates.radial_velocity_correction(kind='barycentric').to(u.km / u.s).value

    return bvc_value


def correct_for_radial_velocity(wave_length, flux, radial_velocity, name_obj,
                                plot=False, merge_type=None):
    """
        Correct wave length array for radial velocity

        Parameters
        ----------
        wave_length         : `numpy.ndarray`
            Wavelength data

        flux                : `numpy.ndarray`
            Flux data

        radial_velocity     : `float`
            Radial velocity

        name_obj            : `string`
            Name of the object

        plot                : `boolean`, optional
            If `True` the plot will be shown
            Default is ``False``.

        merge_type          : `string`, optional
            String characterizing the merging procedure. It will be part of
             the filename.
            Default is ``None``.

        Returns
        -------
        corr_wl             : `numpy.ndarray`
            Redial velocity corrected wavelength
    """
    #   Correct wave length data
    if radial_velocity:
        corrected_wave_length = wave_length / (np.sqrt(
            (c + radial_velocity * 1000) / (c - radial_velocity * 1000)
        ))
    else:
        corrected_wave_length = wave_length

    #   Get plot limits -> use line wavelength + a margin of 20AA
    id_x_min = np.argmin(np.abs(corrected_wave_length - (5889.95 - 20)))
    id_x_max = np.argmin(np.abs(corrected_wave_length - (5895.92 + 20)))

    plt.plot(
        corrected_wave_length[id_x_min:id_x_max],
        flux[id_x_min:id_x_max],
        label="Measured Flux",
        zorder=5,
        color="navy",
    )
    plt.title(f"Vicinity of the sodium 5589 duplett for star {name_obj}")
    plt.axvline(5889.95, color="grey", linestyle="--", zorder=3)
    plt.axvline(5895.92, color="grey", linestyle="--", zorder=3)
    plt.tight_layout()
    if plot:
        plt.show()
    if merge_type is not None:
        plt.savefig(f"output/Vicinity_sodium_{name_obj}_{merge_type}.pdf")
    else:
        plt.savefig(f"output/Vicinity_sodium_{name_obj}.pdf")
    plt.close()

    return corrected_wave_length


def evaluate_flux_difference(scaling_factor, flux_1, flux_2):
    """
        Calculates the number of wavelength points for which the
        flux difference is below a certain threshold, so this can
        be maximized by an external function.

        Parameters
        ----------
        scaling_factor  : `float`
            Initial value for the flux scaling between spectra

        flux_1          : `numpy.ndarray`
            Flux of spectrum 1

        flux_2          : `numpy.ndarray`
            Flux of spectrum 2

        Returns
        -------
                        : `float`
            Inverse of number of flux values below threshold
    """
    #   Set flux threshold
    limit = np.abs(0.02 * np.nanmedian(flux_1))

    #   Calculate points where flux differences is below the threshold
    good_points = (np.abs(flux_1 - flux_2 * scaling_factor) < limit).sum()

    #   Sanitize number of points, since we will return the inverse
    if not good_points:
        good_points = 0.00001

    return 1. / good_points


def normalize_spectrum_interval(wave_length, flux, median_window,
                                normalization_order):
    """
        Normalize a spectrum. If the normalization is not sufficient, the
        spectrum will be split in segments and the normalization will
        be repeated for those segments. The segments will be merged afterward.

        Parameters
        ----------
        wave_length         : `numpy.ndarray`
            Array with wavelength data

        flux                : `numpy.ndarray`
            Array with flux data

        median_window       : `integer`
            Window in Pixel used in median smoothing

        normalization_order : `integer`
            Order of the polynomial that is used for normalization

        Returns
        -------
        wave_length         : `numpy.ndarray`
            Array with normalized wavelength data

        flux                : `numpy.ndarray`
            Array with normalized  flux data
    """
    #   Create Spectrum1D objects
    spectrum_1d = Spectrum1D(spectral_axis=wave_length, flux=flux)

    #   Normalize the merged spectrum
    spectrum_1d, std = normalize_spectrum(spectrum_1d, normalization_order=normalization_order)

    #   Split spectrum in segments, if standard deviation is too high
    if std > 0.05:
        n_segments = 10
        n_wavelength_points = len(wave_length)
        step_size = int(n_wavelength_points / n_segments)
        segments = []

        #   Loop over segments
        i_old = 0
        for i in range(step_size, step_size * n_segments, step_size):
            #   Cut segments and add overlay range to the
            #   segments, so that the normalization afterburner
            #   can take effect
            overlap = int(step_size * 0.15)
            if i == step_size:
                flux_segment = flux[i_old:i + overlap]
                wave_length_segment = wave_length[i_old:i + overlap]
            elif i == n_segments - 1:
                flux_segment = flux[i_old - overlap:]
                wave_length_segment = wave_length[i_old - overlap:]
            else:
                flux_segment = flux[i_old - overlap:i + overlap]
                wave_length_segment = wave_length[i_old - overlap:i + overlap]
            i_old = i

            #   Create Spectrum1D objects for the segments
            segments.append(Spectrum1D(
                spectral_axis=wave_length_segment,
                flux=flux_segment,
            ))

        #   Normalize & merge spectra
        wave_length, flux = normalize_and_merge_spectra(
            segments,
            median_window=median_window,
            normalization_order=normalization_order,
        )

    return wave_length, flux


def normalize_spectrum(spectrum_1d, median_window=61, normalization_order=9):
    """
        Normalize a spectrum

        Parameters
        ----------
        spectrum_1d             : `specutils.Spectrum1D`
            Spectrum to normalize

        median_window           : `int`, optional
            Window in Pixel used in median smoothing
            Default is ``61``.

        normalization_order     : `int`, optional
            Order of the polynomial used to find the continuum
            Default is ``9``.

        Returns
        -------
        normalized_spectrum_1d  : `specutils.Spectrum1D`
            Normalized spectrum
    """
    #   Regions that should not be used for continuum estimation,
    #   such as broad atmospheric absorption bands
    excluded_regions = [
        SpectralRegion(4295. * u.AA, 4315. * u.AA),
        # SpectralRegion(6860.* u.AA, 6880.* u.AA),
        SpectralRegion(6860. * u.AA, 6910. * u.AA),
        # SpectralRegion(7590.* u.AA, 7650.* u.AA),
        SpectralRegion(7590. * u.AA, 7680. * u.AA),
        SpectralRegion(9260. * u.AA, 9420. * u.AA),
        # SpectralRegion(11100.* u.AA, 11450.* u.AA),
        # SpectralRegion(13300.* u.AA, 14500.* u.AA),
        SpectralRegion(4080. * u.AA, 4120. * u.AA),
        SpectralRegion(4320. * u.AA, 4360. * u.AA),
        SpectralRegion(4825. * u.AA, 4885. * u.AA),
        SpectralRegion(6530. * u.AA, 6600. * u.AA),
    ]

    #   First estimate of the continuum
    #   -> will be two for late type stars because of the many absorption lines
    #   -> to limit execution time use simple LinearLSQFitter()
    #       -> reduces normalization accuracy
    _continuum = fit_generic_continuum(
        spectrum_1d,
        model=models.Chebyshev1D(normalization_order),
        fitter=fitting.LinearLSQFitter(),
        median_window=median_window,
        exclude_regions=excluded_regions,
    )(spectrum_1d.spectral_axis)

    #   Normalize spectrum
    normalized_spectrum_1d = spectrum_1d / _continuum

    #   Sigma clip the normalized spectrum to rm spectral lines
    clipped_flux = sigma_clip(
        normalized_spectrum_1d.flux.value,
        sigma_lower=1.25,
        sigma_upper=3.,
        # sigma_upper=4.,
        axis=0,
        grow=1.,
    )

    #   Calculate mask
    mask = np.invert(clipped_flux.recordmask)

    #   Make new spectrum
    masked_spectrum_1d = Spectrum1D(
        spectral_axis=spectrum_1d.spectral_axis[mask],
        flux=spectrum_1d.flux[mask],
    )

    # Determine new continuum
    _continuum = fit_generic_continuum(
        masked_spectrum_1d,
        model=models.Chebyshev1D(normalization_order),
        fitter=fitting.LinearLSQFitter(),
        median_window=median_window,
        exclude_regions=excluded_regions,
    )(spectrum_1d.spectral_axis)

    #   Normalize spectrum again
    normalized_spectrum_1d = spectrum_1d / _continuum

    return normalized_spectrum_1d, mad_std(normalized_spectrum_1d.flux)


def normalize_and_merge_spectra(spectra_list, median_window=61, normalization_order=9):
    """
        Normalize spectra and merge them afterward

        Parameters
        ----------
        spectra_list        : `list` of `specutils.Spectrum1D`
            Spectra to normalize and merge

        median_window       : `int`, optional
            Window in Pixel used in median smoothing
            Default is ``61``.

        normalization_order : `int`, optional
            Polynomial order used for normalization
            Default is ``9``.

        Returns
        -------
                            : `list` of `specutils.Spectrum1D`
            Normalized and merged spectrum
    """
    #   Normalize spectra
    normalized_spectra_list = []
    for spec in spectra_list:
        normalized_spectra_list.append(
            normalize_spectrum(
                spec,
                median_window=median_window,
                normalization_order=normalization_order
            )[0]
        )

    #   Merge spectra
    return merge_normalized_spectra(normalized_spectra_list)


def merge_normalized_spectra(spectra_list):
    """
        Merge normalized spectra

        Parameters
        ----------
        spectra_list            : `list` of `specutils.Spectrum1D`
            List with normalized spectra

        Returns
        -------
        resampled_wavelength    : `list` of `float`
            Common wave length scale

        median_merged_flux      : `list` of `float`
            Merged flux data

    """
    #   Extract wavelength and flux data
    wave_length = []
    flux = []
    for spec in spectra_list:
        wave_length.append(spec.spectral_axis)
        flux.append(spec.flux)

    #   Build common wave length scale
    resampled_wavelength = np.sort(np.concatenate(wave_length))

    #   Prepare list for flux
    resampled_flux = []

    #   Get flux unit
    flux_unit = flux[0].unit

    #   Loop over all spectra
    for i, wave_length_element in enumerate(wave_length):
        #   Remove wavelength duplicates from the input arrays
        _, indices = np.unique(wave_length_element, return_index=True)
        wave_length_element = wave_length_element[indices]
        flux[i] = flux[i][indices]

        #   Interpolate flux on new wavelength grid
        interpolation_object = interpolate.interp1d(
            wave_length_element,
            flux[i],
            kind='cubic',
            bounds_error=False,
        )
        resampled_flux.append(interpolation_object(resampled_wavelength))

    #   Normalization afterburner
    #   -> attempts to improve normalization in the overlap regions
    resampled_flux = normalization_afterburner(resampled_flux)

    #   Merge flux
    median_merged_flux = np.nanmedian(resampled_flux, axis=0)

    return resampled_wavelength, median_merged_flux * flux_unit


def normalization_afterburner(flux_list):
    """
        Afterburner that attempts to correct spectra (especially echelle
        orders) that are not well normalized before merging.
            -> only works if there is an overlap region
            -> assumes that there are only two spectra/orders that overlap
            -> assumes that one of the spectra is well normalized

        Parameters
        ----------
        flux_list       : `list` of `numpy.ndarray`
            Flux of the spectra to normalize

        Returns
        -------
        flux_list       : `list` of `numpy.ndarray`
            Normalized flux
    """

    for i in range(1, len(flux_list)):
        #   Calculate flux difference
        flux_difference = flux_list[i - 1] - flux_list[i]

        #   Check if overlap region exists
        #   -> if not, do nothing
        if ~np.all(np.isnan(flux_difference)):
            #   Determine where flux is negative
            negative_flux_difference_bins = np.argwhere(
                flux_difference < -0.01
            )

            #   Find consecutive elements -> build groups
            negative_flux_difference_ranges = consecutive(
                negative_flux_difference_bins.flatten()
            )

            #   Loop over negative groups
            for negative_difference in negative_flux_difference_ranges:
                #   Restrict to groups with more than 5 elements
                if len(negative_difference) > 5:
                    #   Check if flux of the first spectrum
                    #   is negative in this area
                    #   -> if yes replace
                    #   -> otherwise replace flux in the second spectrum
                    if (np.nanmedian(flux_list[i - 1][negative_difference] - 1)
                            < -0.05):
                        flux_list[i - 1][negative_difference] = flux_list[i][negative_difference]
                    else:
                        flux_list[i][negative_difference] = flux_list[i - 1][negative_difference]

            #   Determine where flux is positive
            positive_flux_difference_bins = np.argwhere(flux_difference > 0.01)

            #   Find consecutive elements -> build groups
            positive_flux_difference_ranges = consecutive(
                positive_flux_difference_bins.flatten()
            )

            #   Loop over positive groups
            for positive_difference in positive_flux_difference_ranges:
                #   Restrict to groups with more than 5 elements
                if len(positive_difference) > 5:
                    #   Check if flux of the second spectrum
                    #   is negative in this area
                    #   -> if yes replace
                    #   -> otherwise replace flux in the first spectrum
                    if (np.nanmedian(flux_list[i][positive_difference] - 1)
                            < -0.05):
                        flux_list[i][positive_difference] = flux_list[i - 1][positive_difference]
                    else:
                        flux_list[i - 1][positive_difference] = flux_list[i][positive_difference]

    return flux_list


def adjust_intersected_spectra(flux_1, flux_2, debug_plot=False, wave_length=None,
                               id_intersection=None):
    """
        Combines two spectra whose flux intersects

        Goal        : Tries to avoid or at least reduce jumps, while merging
        ----

        Parameters
        ----------
        flux_1          : `numpy.ndarray`
            Flux of the first spectrum

        flux_2          : `numpy.ndarray`
            Flux of the second spectrum

        debug_plot      : `boolean`, optional
            If ``True`` the intersection range will be marked on the
            debug plot
            Default is ``False``.

        wave_length     : `numpy.ndarray`, optional
            Wavelength array
            Default is ``None``.

        id_intersection : `integer`, optional
            ID of the current "intersection"
            Default is ``None``.

        Returns
        -------
        flux_1          : `numpy.ndarray`
            Modified flux of the first spectrum

        flux_2          : `numpy.ndarray`
            Modified flux of the second spectrum
    """
    #   Calculate flux difference
    flux_difference = flux_1 - flux_2

    #   Check if overlap exists?
    if ~np.all(np.isnan(flux_difference)):
        #   Calculate range where the flux difference < 10% of maximum of
        #   the flux difference to identify the range where flux intersects
        absolute_flux_difference = np.absolute(flux_difference)
        ids_overlap_region = np.argwhere(
            absolute_flux_difference / np.nanmax(absolute_flux_difference) <= 0.1
        )

        #   Return is no intersection was identified
        if len(ids_overlap_region) == 0:
            return flux_1, flux_2

        #   Find and restrict to center of intersection region
        length_overlap_region = len(ids_overlap_region)
        center_overlap_region = ids_overlap_region[
            int(length_overlap_region / 2.)
        ][0]

        #   Median flux around 10% of intersection region
        #   (plus 1 to ensure that range is not 0)
        intersection_region = int(length_overlap_region * 0.05) + 1
        intersection_region_start = center_overlap_region - intersection_region
        intersection_region_end = center_overlap_region + intersection_region
        median_flux_1_intersection_region = np.median(
            flux_1[intersection_region_start:intersection_region_end]
        )
        if median_flux_1_intersection_region == 0.:
            median_flux_1_intersection_region = 1.

        #   Cut flux range to the overlap range
        flux_difference_overlap = flux_difference[~np.isnan(flux_difference)]

        #   First and last flux point
        start_value_flux_difference_overlap = flux_difference_overlap[0]
        end_value_flux_difference_overlap = flux_difference_overlap[-1]

        #   Find index of the above
        id_start_overlap = np.nanargmin(
            np.absolute(flux_difference - start_value_flux_difference_overlap)
        )
        id_end_overlap = np.nanargmin(
            np.absolute(flux_difference - end_value_flux_difference_overlap)
        )

        #   Plot markers for overlapping region
        if debug_plot:
            plt.plot(
                [
                    wave_length[id_start_overlap].value,
                    wave_length[id_start_overlap].value
                ],
                [
                    0.9 * flux_1[id_start_overlap],
                    1.1 * flux_1[id_start_overlap]
                ],
                color='b',
                linestyle='--',
            )
            plt.plot(
                [
                    wave_length[id_end_overlap].value,
                    wave_length[id_end_overlap].value
                ],
                [
                    0.9 * flux_1[id_end_overlap],
                    1.1 * flux_1[id_end_overlap]
                ],
                color='b',
                linestyle='--',
            )
            plt.plot(
                [
                    wave_length[center_overlap_region].value,
                    wave_length[center_overlap_region].value
                ],
                [
                    0.9 * flux_1[center_overlap_region],
                    1.1 * flux_1[center_overlap_region]
                ],
                color='r',
                linestyle='--',
            )
            plt.text(
                wave_length[center_overlap_region].value,
                1.2 * flux_1[center_overlap_region],
                ' ' + str(id_intersection) + ' ',
                # rotation=90,
                ha='center',
                va='top',
                fontdict=font,
            )

        #   Calculate 3% of the length of the overlap range
        overlap_range_3_percent = int(len(flux_difference_overlap) * 0.03)

        #   Ensure that the above is at least 1
        if overlap_range_3_percent == 0:
            overlap_range_3_percent = 1

        #   Median flux of the first and last 3% in terms of bins
        median_flux_difference_overlap_start = np.median(
            flux_difference_overlap[0:overlap_range_3_percent]
        )
        median_flux_difference_overlap_end = np.median(
            flux_difference_overlap[overlap_range_3_percent * -1:]
        )

        #   If the flux difference is larger than 3% at one edge of the overlap
        #   region, remove the overlapping flux ranges, since the following
        #   merging process would in this case lead to jumps in the merged spectrum
        first_edge_evaluation = np.abs(
            median_flux_difference_overlap_start / median_flux_1_intersection_region
        )
        second_edge_evaluation = np.abs(
            median_flux_difference_overlap_end / median_flux_1_intersection_region
        )
        if first_edge_evaluation > 0.03 or second_edge_evaluation > 0.03:
            flux_2[id_start_overlap:center_overlap_region] = np.nan
            flux_1[center_overlap_region:id_end_overlap] = np.nan

    return flux_1, flux_2


def match_two_spectra(flux_1, flux_2, initial_normalization_factor):
    """
        Normalize and adjust flux of two spectra

        Parameters
        ----------
        flux_1                          : `numpy.ndarray`
            Flux of first spectrum

        flux_2                          : `numpy.ndarray`
            Flux of second spectrum

        initial_normalization_factor    : `float`
            Initial value for the minimization algorithm used to estimate
            the flux difference between `flux_1` and `flux_2`. Usually the
            normalization factor of previous iteration/order

        Returns
        -------
                                        : `numpy.ndarray`
            Flux of second spectrum (`flux_2`) scaled to the first
            one (`flux_1`)

        norm_factor                     : `float`
            Scaling factor of `flux_2`

        Idea: https://stackoverflow.com/questions/13846213/create-composite-spectrum-from-two-unnormalized-spectra
    """

    #   Calculate normalization factor between `flux_1` and `flux_2`, using
    #   scipy.optimize and a normalization function `evaluate_flux_difference`
    optimization_result = optimize.basinhopping(
        evaluate_flux_difference,
        initial_normalization_factor,
        1000,
        minimizer_kwargs={'args': (flux_1, flux_2)}
        # )
    )
    #   If the minimization algorithm fails, use the normalization factor of
    #   previous order
    if optimization_result.success:
        norm_factor = optimization_result.x
    else:
        norm_factor = initial_normalization_factor

    #   Adjust flux and return
    return flux_2 * norm_factor, norm_factor


def merge_spectra(wave_length, flux, trim_value_start_order=400,
                  debug_plot=False):
    """
        Resample spectra on a common wavelength scale and merge those

        Parameters
        ----------
        wave_length             : `list` of `numpy.ndarray`'s
            List of the wavelength ranges

        flux                    : `list` of `numpy.ndarray`'s
            List of the flux ranges corresponding to the wavelength ranges
            in `wave` that should be merged

        trim_value_start_order  : `integer`, optional
            The number of wavelength bins that should be removed from the start
            of the orders.
            Default is ``400``.

        debug_plot              : `boolean`, optional
            If ``True`` the intersection range will be marked on the
            debug plot
            Default is ``False``.

        Returns
        -------
                                : `numpy.ndarray`
            The resampled wavelength scale and the merged flux

                                : `numpy.ndarray`
            The resampled and merged flux
    """
    #   Build common wavelength scale
    resampled_wavelength = np.unique(np.sort(np.concatenate(wave_length)))

    #   Prepare list for flux
    resampled_flux = []

    #   Initial normalization factor
    initial_normalization_factor = 1.

    #   Get flux unit
    flux_unit = flux[0].unit

    print(f"{Bcolors.BOLD}   Match order fluxes ...{Bcolors.ENDC}")
    for i, wave in enumerate(wave_length):
        #   Calculate positions where flux is 0. and rm those
        index_null_flux = np.argwhere(flux[i] == 0.)

        wave_cleaned = np.delete(wave, index_null_flux)
        flux_cleaned = np.delete(flux[i], index_null_flux)

        #   Trim start of the spectral ranges to improve normalization
        #   This proves to be useful in some cases and does not do harm in
        #   the other.
        if i == 0:
            wave_cleaned = wave_cleaned[600:]
            flux_cleaned = flux_cleaned[600:]
        if i > 0:
            wave_cleaned = wave_cleaned[trim_value_start_order:]
            flux_cleaned = flux_cleaned[trim_value_start_order:]

        #   Interpolate flux on common wavelength scale
        interpolation_object = interpolate.interp1d(
            wave_cleaned,
            flux_cleaned,
            kind='cubic',
            bounds_error=False,
        )
        resampled_flux.append(interpolation_object(resampled_wavelength))

        if i > 0:
            print(f"   Adjust flux of order {i} to {i + 1}\r", end="")
            #   Adjust flux of the individual spectra
            resampled_flux[i], initial_normalization_factor = match_two_spectra(
                resampled_flux[i - 1],
                resampled_flux[i],
                initial_normalization_factor,
            )
        if debug_plot:
            plt.step(resampled_wavelength, resampled_flux[i])

    #   Improve overlapping edges before merging
    for i, wave in enumerate(wave_length):
        if i > 0:
            print(
                f"   Improve overlapping edges for order {i} and {i + 1}\r",
                end=""
            )
            #   Manipulate spectra if they "intersect"
            resampled_flux[i - 1], resampled_flux[i] = adjust_intersected_spectra(
                resampled_flux[i - 1],
                resampled_flux[i],
                debug_plot=debug_plot,
                wave_length=resampled_wavelength,
                id_intersection=i,
            )
    print("                                                        \r", end="")

    #   Merge flux and remove residuals NANS
    print(f"{Bcolors.BOLD}   Merge orders... {Bcolors.ENDC}")
    median_merged_flux = np.nanmedian(resampled_flux, axis=0)
    index_nan_flux = np.argwhere(np.isnan(median_merged_flux))

    return (np.delete(resampled_wavelength, index_nan_flux),
            np.delete(median_merged_flux, index_nan_flux) * flux_unit)


def normalize_spectrum_fabian(wave_length, flux, object_name="",
                              lines_to_skip=None,
                              range_to_skip_around_lines_to_skip=None):
    """
        Normalize a spectrum to its continuum

        Parameters
        ----------
        wave_length                         : `numpy.ndarray`
            Wavelength data

        flux                                : `numpy.ndarray`
            Flux data

        object_name                         : `string`, optional
            Name of the object
            Default is an empty string.

        lines_to_skip                       : `list` or None, optional
            Line wavelength that should be skipped during continuum
            determination
            Default is ``None``.

        range_to_skip_around_lines_to_skip  : `list` or `None`, optional
            Wavelength limits that should be applied around each line that
            should be skipped
            Default is ``None``.

        Returns
        -------
        normalized_flux                     : `numpy.ndarray`
            Normalized flux data

        continuum                           : `numpy.ndarray`
            Continuum data points
    """
    #   Sanitize variables
    if lines_to_skip is None:
        lines_to_skip = []
    if range_to_skip_around_lines_to_skip is None:
        range_to_skip_around_lines_to_skip = []

    #   Calculate median of the step size between wavelength points
    wave_length_step = np.median(np.diff(wave_length))

    medium_window_size = 0.5
    max_window_size = 12
    #   Calculate number of elements in the maximum and median window
    elements_in_medium_window = math.floor(
        medium_window_size / wave_length_step
    )
    elements_in_max_window = math.floor(max_window_size / wave_length_step)

    if elements_in_max_window == 0 or elements_in_medium_window == 0:
        raise ValueError(
            "Medium/Maximum window sizes need to be bigger than a "
            "wavelength step!"
        )

    #   Apply median and maximum filter with the corresponding window sizes
    flux_for_interpolation = median_filter(
        flux,
        size=elements_in_medium_window,
    )
    flux_for_interpolation = maximum_filter(
        flux_for_interpolation,
        size=elements_in_max_window,
    )

    #   Cut spectral lines that would spoil the continuum fit such as broad
    #   hydrogen lines
    spectral_lines_mask = np.ones(np.shape(wave_length), dtype=bool)
    for i, line_wave_length in enumerate(lines_to_skip):
        line_mask = np.logical_or(
            wave_length > line_wave_length + range_to_skip_around_lines_to_skip[i],
            wave_length < line_wave_length - range_to_skip_around_lines_to_skip[i]
        )
        spectral_lines_mask = np.logical_and(
            spectral_lines_mask,
            line_mask
        )

    flux_for_interpolation = flux_for_interpolation[spectral_lines_mask]
    wave_length_for_interpol = wave_length[spectral_lines_mask]

    #   Interpolate to find the continuum
    continuum_fit = InterpolatedUnivariateSpline(
        wave_length_for_interpol,
        flux_for_interpolation,
        k=2,
    )

    #   Normalize spectrum
    continuum = continuum_fit(wave_length)
    normalized_flux = flux / continuum

    if object_name != "":
        np.savetxt(
            f"flux_normalized_{object_name}.csv",
            np.transpose((wave_length, normalized_flux)),
            delimiter=",",
        )

    return normalized_flux, continuum


def read_baches(file_name):
    """
        Read baches FITS files with individual orders created by MIDAS

        Parameters
        ----------
        file_name           : `string`
            Name of the file to read

        Returns
        -------
        wave_length_list    : `list of `numpy.ndarray`
            List with the wavelength data of the individual orders

        flux_list           : `list of `numpy.ndarray`
            List with the flux data of the individual orders
    """
    file = fits.open(file_name)
    file_data = file[0].data

    file_header = file[0].header
    history = file_header['HISTORY']
    # ref_pixel        = FileHeader['CRPIX1']
    # coord_ref_pixel  = FileHeader['CRVAL1']
    wave_per_pixel = file_header['CDELT1']

    #   Extract start points from HISTORY section of the FITS Header
    k = 1
    start_points = []
    for line in history:
        line_list = line.split()
        if not line_list:
            k = 1
            continue
        if k == 0:
            start_points = start_points + line_list
        if k == 1:
            if 'WSTART' not in line_list[0]:
                continue
            else:
                k = 0
    start_points = np.array(start_points, dtype=float)

    #   Convert fluxes to numpy array and calculate wavelength points for
    #   all flux values
    flux_array = np.array(file_data, dtype=float)
    wave_length = (np.ones(flux_array.shape) *
                   np.arange(0, flux_array.shape[1]) * wave_per_pixel)
    wave_length = wave_length + start_points.reshape(start_points.size, 1)

    #   Limit to flux values != 0
    wave_length_list = []
    flux_list = []
    for i, flux in enumerate(flux_array):
        mask_zero_flux = np.nonzero(flux)
        wave_length_list.append(wave_length[i, mask_zero_flux][0] * u.AA)
        flux_list.append(flux[mask_zero_flux] * u.adu)

    return wave_length_list, flux_list


def read_baches_merged(file_name):
    """
        Read baches FITS file with spectrum merged by MIDAS

        Parameters
        ----------
        file_name           : `string`
            Name of the file to read

        Returns
        -------
        wave_length_merged  : `numpy.ndarray`
            Wavelength data of the spectrum merged by MIDAS

        flux_merged         : `numpy.ndarray`
            Flux data of the spectrum merged by MIDAS
    """
    file_merged = fits.open(file_name)
    file_data_merged = file_merged[0].data

    file_header_merged = file_merged[0].header
    # ref_pixel_merged       = FileHeader_merged['CRPIX1']
    reference_pixel_merged = file_header_merged['CRVAL1']
    wave_per_pixel_merged = file_header_merged['CDELT1']

    #   Number of lines in the file
    file_data_lines_merged = len(file_data_merged)

    wave_length_merged = []
    for i in range(0, file_data_lines_merged):
        wave_length_merged.append(
            reference_pixel_merged + i * wave_per_pixel_merged
        )

    start_cut = 200
    end_cut = 200
    wave_length_merged = wave_length_merged[start_cut:-end_cut]
    flux_merged = file_data_merged[start_cut:-end_cut]

    return wave_length_merged * u.AA, flux_merged * u.adu


def check_file_name(file_name):
    """
        Check that the file has the correct ending

        Parameters
        ----------
        file_name            : `string`
            Path to the file

        Returns
        -------
        file_name            : `string`
            Modified path to the file
    """
    if os.path.isfile(file_name) is False:
        if file_name.find(".fit") > 0:
            file_name = file_name.replace(".fit", ".FIT")
            print(f"change .fit to .FIT in {file_name}")
        elif ".FIT" in file_name:
            file_name = file_name.replace(".FIT", ".fit")
            print(f"change .FIT to .fit in {file_name}")
        else:
            print(
                f"{Bcolors.FAIL}!!!!! Check the right spelling of the file "
                f"extension and the filepath !!!!!{Bcolors.ENDC}"
            )

    return file_name


def add_idents(wave_length, flux, line_file, min_flux=None, max_flux=None):
    """
        Add idents to the plot

        Parameters
        ----------
        wave_length      : `numpy.ndarray`
            Wavelength array

        flux             : `numpy.ndarray`
            Flux array

        line_file        : `string`
            Name of the ident file

        min_flux         : `float`, optional
            Minimum plot range on the Y axis
            Default is ``None``.

        max_flux         : `float`, optional
            Maximum plot range on the Y axis
            Default is ``None``.
    """
    #   Defining plot range
    if min_flux is None:
        min_flux = min(flux)
        min_flux = max(0, min_flux)
    if max_flux is None:
        max_flux = max(flux)

    flux_offset = (max_flux - min_flux) * 0.01
    # x_offset = (max(wave_length) - min(wave_length)) * 0.01

    plot_minimum_y = min_flux - flux_offset
    plot_maximum_y = max_flux + flux_offset

    #   Setting plot positions for ident lines
    plot_height = plot_maximum_y - plot_minimum_y
    plot_middle_y = (plot_minimum_y + plot_maximum_y) / 2.0
    plot_upper_y = plot_minimum_y + 0.90 * plot_height
    plot_lower_y = plot_minimum_y + 0.10 * plot_height
    plot_upper_y_cut_1 = plot_minimum_y + 0.80 * plot_height
    plot_lower_y_cut_1 = plot_minimum_y + 0.20 * plot_height
    plot_upper_y_cut_2 = plot_minimum_y + 0.78 * plot_height
    plot_lower_y_cut_2 = plot_minimum_y + 0.22 * plot_height

    #   Interpolate on data to find point for ident
    flux_interpolation_fit = interpolate.interp1d(wave_length, flux)

    #   Open ident file
    try:
        lines = open(line_file, "r")
    except RuntimeError as e:
        print(
            f"{Bcolors.FAIL}   Line file not found. Check variable "
            f"'lineFile'. Specified was {line_file}. {Bcolors.ENDC}"
        )
        raise e

    #   Plot idents to figure
    for line in lines:
        line_list = line.split()

        if len(line_list) == 1:
            print(
                f"{Bcolors.WARNING}   [WARNING] Broken identification found "
                f"as '{line}', must consist of [wavelength(s) + name]. I will "
                f"skip this one.{Bcolors.ENDC}"
            )
            continue
        try:
            float(line_list[0])
        except ValueError:
            print(
                f"{Bcolors.WARNING}   [WARNING] Broken identification found "
                f"as '{line}', first entry not a number. I will skip this "
                f"one.{Bcolors.ENDC}"
            )
            continue

        #   Single ident
        if len(line_list) == 2:
            ident_line = line_list

            #   Only plot if in range
            if min(wave_length) <= float(ident_line[0]) <= max(wave_length):
                #   Plot ident point in upper or lower half of figure
                if flux_interpolation_fit(ident_line[0]) >= plot_middle_y:
                    plt.plot(
                        [float(ident_line[0]), float(ident_line[0])],
                        [
                            flux_interpolation_fit(float(ident_line[0])),
                            plot_lower_y_cut_2
                        ],
                        color='r',
                        linestyle='-',
                    )
                    plt.text(
                        float(ident_line[0]),
                        plot_lower_y_cut_2,
                        ' ' + ident_line[1] + ' ',
                        rotation=90,
                        ha='center',
                        va='top',
                        fontdict=font,
                    )
                else:
                    plt.plot(
                        [float(ident_line[0]), float(ident_line[0])],
                        [
                            plot_upper_y_cut_2,
                            flux_interpolation_fit(float(ident_line[0]))
                        ],
                        color='r',
                        linestyle='-',
                    )
                    plt.text(
                        float(ident_line[0]),
                        plot_upper_y_cut_2,
                        ' ' + ident_line[1] + ' ',
                        rotation=90,
                        ha='center',
                        va='bottom',
                        fontdict=font,
                    )

        #   Multi ident
        if len(line_list) > 2:
            wave_length_line_elements = []
            flux_line_elements = []
            for wave_length_line_element in line_list[:-1]:
                wave_length_line_elements.append(
                    float(wave_length_line_element)
                )
            wave_length_line_center = (sum(wave_length_line_elements) /
                                       float(len(wave_length_line_elements))
                                       )
            ident_name = str(line_list[-1:])
            if (max(wave_length_line_elements) <= max(wave_length) and
                    min(wave_length_line_elements) >= min(wave_length)):
                for wave_length_line_element in line_list[:-1]:
                    flux_line_elements.append(
                        flux_interpolation_fit(float(wave_length_line_element))
                    )
                flux_center = (sum(flux_line_elements) /
                               float(len(wave_length_line_elements))
                               )

                if flux_center <= plot_middle_y:
                    plt.plot(
                        [wave_length_line_center, wave_length_line_center],
                        [plot_upper_y, plot_upper_y_cut_1],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        wave_length_line_center,
                        plot_upper_y,
                        ' ' + ident_name[2:-2] + ' ',
                        rotation=90,
                        ha='center',
                        va='bottom',
                        fontdict=font,
                    )
                    for element in wave_length_line_elements:
                        plt.plot(
                            [element, element],
                            [
                                flux_interpolation_fit(element), 
                                plot_upper_y_cut_2
                            ],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )
                        plt.plot(
                            [wave_length_line_center, element],
                            [plot_upper_y_cut_1, plot_upper_y_cut_2],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )

                if flux_center > plot_middle_y:
                    plt.plot(
                        [
                            wave_length_line_center, 
                            wave_length_line_center
                        ],
                        [plot_lower_y, plot_lower_y_cut_1],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        wave_length_line_center,
                        plot_lower_y,
                        ' ' + ident_name[2:-2] + ' ',
                        rotation=90,
                        ha='center',
                        va='top',
                        fontdict=font,
                    )
                    for element in wave_length_line_elements:
                        plt.plot(
                            [element, element],
                            [
                                flux_interpolation_fit(element), 
                                plot_lower_y_cut_2
                            ],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )
                        plt.plot(
                            [wave_length_line_center, element],
                            [plot_lower_y_cut_1, plot_lower_y_cut_2],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )
    lines.close()


def add_idents_fabian(lines_dict, x_plotting_limits, current_subplot):
    """
        Add idents to the plot

        Parameters
        ----------
        lines_dict                   : `dictionary`
            Line information: key=line identifier, value[0]=wavelength,
                                                   value[1]=alignment parameter

        x_plotting_limits            : `tuple`
            Limits for the plot in X direction

        current_subplot              : `matplotlib.pyplot.subplots`
            Plot to which the idents should be added.
    """
    for line_identifier, line_locations in lines_dict.items():
        line_wave_lengths = [location[0] for location in line_locations]
        #   Restrict to line within plot range
        mask = np.logical_and(
            line_wave_lengths > x_plotting_limits[0],
            line_wave_lengths < x_plotting_limits[1]
        )
        if sum(mask) != 0:
            for location in line_locations:
                #   Align parameter
                align_parameter = location[1]
                #   Wave length
                line_wave_length = location[0]

                if not x_plotting_limits[0] < line_wave_length < x_plotting_limits[1]:
                    continue

                #   Plot identifier
                x_axis_transform = current_subplot.get_xaxis_transform()
                current_subplot.annotate(
                    line_identifier,
                    xy=(line_wave_length, 1.05),
                    xycoords=x_axis_transform,
                    ha=align_parameter,
                )
                current_subplot.axvline(
                    line_wave_length,
                    color="lightgrey",
                    linestyle="--",
                    zorder=1,
                )


def plot_merged(wave_length_merged, flux_merged, object_name, normalize,
                merge_type=''):
    """
        Plot merged spectrum

        Parameters
        ----------
        wave_length_merged      : `numpy.ndarray`
            Wavelength data

        flux_merged             : `numpy.ndarray`
            Flux data

        object_name             : `string`
            Name of the object

        normalize               : `boolean`
            If `True`, it is assumed that the flux is normalized.

        merge_type              : `string`, optional
            String characterizing the merging procedure. It will be part of
             the filename.
            Default is ``''``.

    """
    print(
        f'      Plot total spectrum:{Bcolors.OKBLUE} {object_name}'
        f'{Bcolors.ENDC}'
    )

    #   Define figure
    figure = plt.figure(figsize=(10, 5))

    #   Set plot range
    flux_offset = (max(flux_merged) - min(flux_merged)) * 0.02
    plt.ylim([min(flux_merged) - flux_offset, max(flux_merged) + flux_offset])
    wave_length_offset = (max(wave_length_merged) - min(wave_length_merged)) * 0.02
    plt.xlim([
        min(wave_length_merged) - wave_length_offset,
        max(wave_length_merged) + wave_length_offset
    ])

    #   Set title and label etc.
    plt.suptitle(
        f"{object_name.replace('_', ' ')} (orders merged by {merge_type})"
    )
    plt.xlabel(r'Wavelength $\lambda\,[\AA]$')
    if normalize:
        plt.ylabel('Normalized flux')

        #   Set plot range
        min_flux = np.max([np.min(flux_merged) * 0.94, 0])
        max_flux = np.min([np.max(flux_merged) * 1.06, 2])
        plt.ylim([min_flux, max_flux])
    else:
        plt.ylabel('Relative flux')
        min_flux = np.max([np.min(flux_merged) * 0.94, 0])
        max_flux = np.max(flux_merged) * 1.06
        plt.ylim([min_flux, max_flux])

    plt.tick_params(top=True, right=True, which='both', direction='in')
    plt.minorticks_on()

    #   Plot spectrum
    plt.plot(wave_length_merged, flux_merged, 'b-', linewidth=0.5)

    #   Save spectrum as PDF file
    print(
        f'      Create spectrum plot{Bcolors.OKBLUE} output/'
        f'spectrum_total_{merge_type}-merged_{object_name}.pdf {Bcolors.ENDC}'
    )
    os.system("mkdir -p output")
    pp = PdfPages(
        f'output/spectrum_total_{merge_type}-merged_{object_name}.pdf'
    )

    pp.savefig(figure, dpi=300, transparent=True)
    pp.close()
    plt.clf()
    plt.close()

    #   Write data to CSV file
    print(
        f'      Write data to{Bcolors.OKBLUE} output/{object_name}_'
        f'{merge_type}-merged_spectrum_total.csv {Bcolors.ENDC}'
    )

    np.savetxt(
        f"output/{object_name}_{merge_type}-merged_spectrum_total.csv",
        np.transpose((wave_length_merged, flux_merged)),
        delimiter=",",
    )
    np.savetxt(
        f"output/{object_name}_{merge_type}-merged_spectrum_total.dat",
        np.transpose((wave_length_merged, flux_merged)),
        delimiter=" ",
    )


def to_roman(number):
    """
        Convert Arabic to Roman numbers

        Parameters
        ----------
        number          : `string`
            Arabic number to convert

        Returns
        -------
                        : `string`
            Roman number
    """
    roman = {
        "1": "I",
        "2": "II",
        "3": "III",
        "4": "IV",
        "5": "V",
        "6": "VI",
        "7": "VII",
        "8": "VIII",
        "9": "IX",
        "10": "X",
        "11": "XI",
        "12": "XII",
        "13": "XIII",
        "14": "XIV",
        "15": "XV"
    }
    return roman[number]


def replace_arabic_with_roman_numbers(string):
    """
        Replaces the Arabic numbers in a string (second part, delimiter=' ')
        with Roman numbers.

        Parameters
        ----------
        string          : `string`
            String to split.

        Returns
        -------
        string          : `string`
            String with Roman numbers.
    """
    if " " in string.strip():
        string_split = string.split(" ")
        return f'{string_split[0]}{to_roman(string_split[1])}'.strip()
    else:
        return string


def generate_lines_from_file(wave_length_array, flux, ion_list,
                             line_file="atomic_lines.tsv",
                             percentage_line_flux_must_be_below_continuum=3.):
    """
        Read line file and prepare line_dict of the ions that are requested.
        Gets the line_dict with the larges gaunt factors. The line_dict need to be
        within the plot range.

        Parameters
        ----------
        wave_length_array                               : `numpy.ndarray`
            Wavelength data

        flux                                            : `numpy.ndarray`
            Flux data

        ion_list                                        : `list` of `string`
            Ions whose line_dict are to be extracted from the file.

        line_file                                       : `string`, optional
            Name/Path to the file with spectral line data
            Default is ``atomic_lines.tsv``.

        percentage_line_flux_must_be_below_continuum    : `float`, optional
            Percent the line flux must be lower than the continuum
            Default is ``3``.

        Returns
        -------
        line_dict                                       : `dict`
            Line information: key=ion identifier, value=wavelength
    """
    #   Interpolate on wavelength and flux to allow determination of flux
    #   level at the individual line positions
    flux_interpolation_fit = interp1d(
        wave_length_array,
        flux,
        fill_value=1,
        bounds_error=False,
    )

    #   Read line file and restrict data to a reasonable range
    try:
        line_data = pd.read_csv(line_file, delimiter="\t")
    except:
        print(
            f"{Bcolors.WARNING}   Line file not found. Check variable "
            f"'lineFile'. Specified was {line_file}. {Bcolors.ENDC}"
        )
        return {}
    line_data = line_data.loc[
        (line_data["wave_A"] > 4200) & (line_data["wave_A"] < 7600)
        ]

    #   Replace Arabic numbers with Roman numbers
    line_data['element'] = line_data['element'].apply(
        replace_arabic_with_roman_numbers
    )

    #   Setup lists and a dict for the line data
    line_dict = {}
    line_wave_length_list = []
    line_strength_list = []
    line_ion_list = []

    #   Restrict line data to strong line_dict and line_dict with the largest
    #   gaunt factors
    for ion in ion_list:
        #   Restrict line data to current ion and the 25 strongest lines
        line_data_subset = line_data.loc[line_data["element"] == ion]
        line_data_subset = line_data_subset.sort_values('loggf').tail(25)

        for _, row in line_data_subset.iterrows():
            #   Get line wavelength and index of closest wavelength in wave
            #   length array
            line_wave_length = float(row["wave_A"])
            line_index = np.argmin(
                np.abs(wave_length_array - line_wave_length)
            )

            #   Mean flux around line wavelength
            try:
                mean_flux_around_line = np.mean(
                    flux[line_index - 40:line_index + 40]
                )
            except IndexError:
                continue

            #   Skip weak (not deep) lines with except for of hydrogen lines
            line_flux_criterion = ((100 - percentage_line_flux_must_be_below_continuum)
                                   / 100 * mean_flux_around_line)
            if flux_interpolation_fit(line_wave_length) > line_flux_criterion:
                if row["element"] != 'HI':
                    continue

            #   Remove/skip lines with the same wavelength if there is one
            #   with a larger gaunt factor
            if line_wave_length in line_wave_length_list:
                line_id = line_wave_length_list.index(line_wave_length)
                if row["loggf"] > line_strength_list[line_id]:
                    line_dict[line_ion_list[line_id]].remove(
                        [line_wave_length, "center"]
                    )

                    line_wave_length_list.pop(line_id)
                    line_strength_list.pop(line_id)
                    line_ion_list.pop(line_id)
                else:
                    continue

            #   Check if ion already exists in the `line_dict` dictionary
            try:
                line_dict[row["element"]]
            except KeyError:
                line_dict[row["element"]] = []

            #   Fill list with line data
            line_dict[row["element"]].append([line_wave_length, "center"])
            line_wave_length_list.append(line_wave_length)
            line_strength_list.append(float(row["loggf"]))
            line_ion_list.append(row["element"])

    #   Convert lists to arrays
    line_wave_length_array = np.array(line_wave_length_list)
    line_strength_array = np.array(line_strength_list)
    ion_array = np.array(line_ion_list)

    #   Remove line_dict that are too close to each other - keep strongest
    for line_wave_length in line_wave_length_array:
        #   Find close lines
        mask = np.logical_and(
            line_wave_length - 4 < line_wave_length_array,
            line_wave_length_array < line_wave_length + 4
        )
        if sum(mask) > 1:
            #   Hydrogen is always the strongest line
            line_index = np.argwhere(line_wave_length_array == line_wave_length)
            if "HI" in ion_array[line_index]:
                strongest_line = line_wave_length
            else:
                #   Find strong line based on gaunt factor
                index_strongest_line = np.argmax(line_strength_array[mask])
                strongest_line = line_wave_length_array[mask][index_strongest_line]

            #   Get other line_dict (weaker line_dict)
            other_lines_wave_length = line_wave_length_array[mask]
            other_lines_wave_length = other_lines_wave_length[other_lines_wave_length != strongest_line]

            #   Remove other line_dict
            for other_line_wave_length in other_lines_wave_length:
                index_other_line = np.argwhere(
                    line_wave_length_array == other_line_wave_length
                )[0][0]
                line_dict[ion_array[index_other_line]].remove(
                    [line_wave_length_array[index_other_line], "center"]
                )

                line_wave_length_array = np.delete(
                    line_wave_length_array,
                    index_other_line
                )
                # line_strength_list = np.delete(line_strength_array, index_other_line)
                ion_array = np.delete(ion_array, index_other_line)

    #   Set horizontal alignment for identifier string, if lines are too close
    for ion, line_list in line_dict.items():
        for line in line_list:
            #   Line wave length
            line_wave_length = line[0]

            #   Find close line_dict
            # if sum(np.logical_and(line_wave_length_array > l - 0.5, line_wave_length_array < l + 0.5)) > 1:
            # np.delete(line_wave_length_array, line_wave_length_array == l)
            # line_dict[element].remove([l, "center"])
            # continue

            #   Identify close line_dict
            mask = np.logical_and(
                line_wave_length_array > line_wave_length - 10,
                line_wave_length_array < line_wave_length + 10
            )
            if sum(mask) > 1:
                other_lines_wave_length = line_wave_length_array[mask]
                other_lines_wave_length = other_lines_wave_length[
                    other_lines_wave_length != line_wave_length
                ]

                #   Set alignment according to line wavelength
                for other_line_wave_length in other_lines_wave_length:
                    line_index = np.argwhere(
                        line_wave_length_array == other_line_wave_length
                    )[0][0]
                    other_ion = ion_array[line_index]
                    other_ion_wave_lengths = np.array(line_dict[other_ion])[:, 0]
                    index_wave_length = list(other_ion_wave_lengths).index(
                        str(other_line_wave_length)
                    )
                    if other_line_wave_length > line_wave_length:
                        line_dict[other_ion][index_wave_length][1] = "left"
                    else:
                        line_dict[other_ion][index_wave_length][1] = "right"

    return line_dict


def plot_panels(wave_length_merged, flux_merged, object_name, normalize,
                line_file, panel_wave_length_range, merge_type=''):
    """
        Plot merged data in individual panels and create PDFs

        Parameters
        ----------
        wave_length_merged             : `numpy.ndarray`
            Wavelength data

        flux_merged             : `numpy.ndarray`
            Flux data

        object_name                : `string`
            Name of the object

        normalize               : `boolean`
            If `True`, it is assumed that the flux is normalized.

        line_file               : `string`
            Path to the file with the line identifications

        panel_wave_length_range        : `float` or `integer`
            Wavelength range for the individual panels

        merge_type              : `string`, optional
            String characterizing the merging procedure. It will be part of
             the filename.
            Default is ``''``.
    """
    print(
        f"      Plot individual panels:{Bcolors.OKBLUE} {object_name}"
        f"{Bcolors.ENDC}"
    )

    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    n_wave_length_steps = len(wave_length_merged)
    i = 1
    start_wave_length_range = 0
    for j in range(0, n_wave_length_steps):
        if (wave_length_merged[j] >= wave_length_merged[0] + panel_wave_length_range * i
                or j + 1 == n_wave_length_steps):

            #   Set wave and flux range
            wave_length_range = wave_length_merged[start_wave_length_range:j]
            flux_range = flux_merged[start_wave_length_range:j]

            #   Define figure and set labels etc.
            plt.figure(figsize=(12, 6))
            plt.suptitle(
                f"{object_name.replace('_', ' ')} ({int(wave_length_range[0])}"
                fr" - {int(wave_length_range[-1])}$\,\AA$, orders merged by {merge_type})"
            )
            plt.xlabel(r'Wavelength $\lambda\,[\AA]$')
            if normalize:
                plt.ylabel('Normalized flux')

                #   Set plot range
                min_flux = np.max([np.min(flux_range) * 0.94, 0])
                max_flux = np.min([np.max(flux_range) * 1.06, 2])
                plt.ylim([min_flux, max_flux])
            else:
                plt.ylabel('Relative flux')
                min_flux = np.max([np.min(flux_range) * 0.94, 0])
                max_flux = np.max(flux_range) * 1.06
                plt.ylim([min_flux, max_flux])

            plt.tick_params(top=True, right=True, which='both', direction='in')
            plt.minorticks_on()

            #   Plot idents
            add_idents(
                wave_length_range, 
                flux_range, 
                line_file, 
                min_flux, 
                max_flux,
            )

            #   Plot spectrum
            plt.step(wave_length_range, flux_range, 'b-', linewidth=0.5)

            #   Save panel as PDF file
            plt.savefig(
                f'{temp_dir.name}/{object_name}_{merge_type}-merged_{i}.pdf',
                bbox_inches='tight',
            )
            plt.close()

            #   Write data to CSV file
            print(
                f'      Write data to{Bcolors.OKBLUE} output/{object_name}_'
                f'{merge_type}-merged_spectrum_{i}.csv {Bcolors.ENDC}'
            )
            np.savetxt(
                f"output/{object_name}_{merge_type}-merged_spectrum_{i}.csv",
                np.transpose((wave_length_range, flux_range)),
                delimiter=",",
            )

            i += 1
            start_wave_length_range = copy.deepcopy(j)

    #   Merge individual plotted bins with pdfunite
    print(
        f"      Create{Bcolors.OKBLUE} output/spectrum_panels_{merge_type}"
        f"-merged_{object_name}.pdf {Bcolors.ENDC}"
    )

    merger = PdfMerger()
    file_path = f"output/spectrum_panels_{merge_type}-merged_{object_name}.pdf"
    for file_ in os.listdir(temp_dir.name):
        merger.append(f"{temp_dir.name}/{file_}")
    with open(file_path, "wb") as new_file:
        merger.write(new_file)

    merger.close()


def plot_panels_fabian(wave_length, flux, object_name, normalize, 
                       lines_to_mark, n_panels=21, merge_type='', 
                       radial_velocity=0., continuum=None, 
                       n_panels_per_page=5):
    """
        Plot merged data in individual panels and create PDFs

        Parameters
        ----------
        wave_length             : `numpy.ndarray`
            Wavelength data

        flux                    : `numpy.ndarray`
            Flux data

        object_name             : `string`
            Name of the object

        normalize               : `boolean`
            If `True`, it is assumed that the flux is normalized.

        lines_to_mark           : `dictionary`
            Line information: key=line identifier, value[0]=wavelength,
                                                   value[1]=alignment parameter

        n_panels                : `integer`, optional
            Number of panels that should be used to display the spectrum
            Default is ``21``.

        merge_type              : `string`, optional
            String characterizing the merging procedure. It will be part of
             the filename.
            Default is ``''``.

        radial_velocity         : `float`, optional
            Radial velocity

        continuum               : `numpy.ndarray`
            Flux data of the continuum

        n_panels_per_page       : `integer`, optional
            Number of panels on each page/plot.
            Default is ``5``.
    """
    print(
        f"      Plot individual panels:{Bcolors.OKBLUE} {object_name}"
        f"{Bcolors.ENDC}"
    )

    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Limit to a range which can be divided by ``n_panels``
    while len(flux) % n_panels != 0:
        wave_length = wave_length[:-1]
        flux = flux[:-1]
        if continuum is not None:
            continuum = continuum[:-1]

    #   Split data into ``n_panels`` ranges
    wave_split = np.split(wave_length, n_panels)
    flux_split = np.split(flux, n_panels)
    if continuum is not None:
        cont_split = np.split(continuum, n_panels)
    else:
        cont_split = [0 for _ in wave_split]

    #   Initialize plot
    figure, axes = plt.subplots(n_panels_per_page, 1, figsize=(8.2, 11.6))

    for i, [wave_length_element, flux_element, continuum_element] in (
            enumerate(zip(wave_split, flux_split, cont_split))):
        #   Axis index
        axis_index = i
        while axis_index >= 0:
            axis_index -= n_panels_per_page
        i += 1

        # plt.tick_params(top=True, right=True, which='both', direction='in')
        # plt.minorticks_on()

        #   Plot data
        axes[axis_index].step(
            wave_length_element,
            flux_element,
            color='#0066ff',
            zorder=5,
        )

        #   Plot continuum
        if continuum is not None:
            axes[axis_index].plot(
                wave_length_element,
                continuum_element,
                color='red',
                zorder=5,
            )

        #   Set ranges
        x_limit = (np.amin(wave_length_element), np.amax(wave_length_element))
        axes[axis_index].set_xlim(x_limit)
        if normalize:
            # min_flux = np.max([np.amin(flux_element) * 0.94, 0])
            min_flux = 0.
            max_flux = np.min([np.amax(flux_element) * 1.06, 2])
            axes[axis_index].set_ylim([min_flux, max_flux])
            axes[axis_index].set_ylabel('Normalized flux')
            axes[axis_index].legend(["Measured"], loc="lower right")
            axes[axis_index].axline(
                (x_limit[0], 1.0),
                xy2=(x_limit[1], 1.0),
                color="darkgrey",
                linestyle="-",
                linewidth=0.5,
                zorder=1,
            )
        else:
            min_flux = np.max([np.amin(flux_element) * 0.94, 0])
            max_flux = np.amax(flux_element) * 1.06
            axes[axis_index].set_ylim([min_flux, max_flux])
            axes[axis_index].set_ylabel('Relative flux')
            if continuum is not None:
                axes[axis_index].legend(
                    ["Measured", "Continuum"],
                    loc="lower right"
                )
            else:
                axes[axis_index].legend(["Measured"], loc="lower right")

        #   Set plot title, labels, etc.
        figure.suptitle(
            f"Spectrum of {object_name} ({merge_type}) - Page "
            f"[{math.ceil(i / n_panels_per_page)}/"
            f"{math.ceil(n_panels / n_panels_per_page)}]",
            fontsize=15,
        )

        axes[n_panels_per_page - 1].set_xlabel(r"$\lambda$ [Å]")

        if ((i - 1) % n_panels_per_page) == 0:
            axes[axis_index].annotate(
                f"Corrected for radial velocity "
                fr"$v_{{\mathrm{{rad}}}}={round(radial_velocity, 2)}$ km/s",
                xy=(axes[axis_index].get_xlim()[0], 1.125),
                xycoords=axes[axis_index].get_xaxis_transform(),
                ha="left",
                fontsize=10,
            )

        #   Plot line identifier
        add_idents_fabian(lines_to_mark, x_limit, axes[axis_index])

        #   Save plot
        if i % n_panels_per_page == 0:
            plt.tight_layout()
            plt.savefig(f"{temp_dir.name}/{int(i / n_panels_per_page)}.pdf")
            plt.close()
            figure, axes = plt.subplots(
                n_panels_per_page,
                1,
                figsize=(8.2, 11.6),
            )

            # #   Write data to CSV file
            # print(
            # '         Write data to{} output/{}_{}-merged_spectrum'
            # '_{}.csv {}'.format(
            # bcolors.OKBLUE,
            # objectname,
            # mtype,
            # i,
            # bcolors.ENDC,
            # )
            # )
            # np.savetxt(
            # f"output/{objectname}_{mtype}-merged_spectrum_{i}.csv",
            # np.transpose((wave_range, flux_range)),
            # delimiter=",",
            # )
    plt.close()

    #   Merge individual plotted bins with pdfunite
    print(
        f"      Create{Bcolors.OKBLUE} output/spectrum_panels_{merge_type}-"
        f"merged_{object_name}.pdf {Bcolors.ENDC}"
    )
    merger = PdfMerger()
    file_path = f"output/spectrum_panels_{merge_type}-merged_{object_name}.pdf"
    for flux_element in os.listdir(temp_dir.name):
        merger.append(f"{temp_dir.name}/{flux_element}")
    with open(file_path, "wb") as new_file:
        merger.write(new_file)

    merger.close()


############################################################################
#                                  Main                                    #
############################################################################

if __name__ == '__main__':
    ###
    #   Prepare stuff
    #
    #  Sanitize object name
    object_name = object_name.replace(' ', '_')

    #  Create output directory
    os.system("mkdir -p output")

    ########################################################################
    #                         MIDAS merged spectrum                        #
    ########################################################################
    print(
        f"{Bcolors.BOLD}Processing spectrum merged by MIDAS{Bcolors.ENDC}"
    )

    ###
    #   Load spectra
    #
    #   Select FIT or fit
    file_with_merged_spectrum = check_file_name(file_with_merged_spectrum)

    #   Open fits file and read the data section and header section into a
    #   2-dimensional array
    print(f"{Bcolors.BOLD}   Read files ...{Bcolors.ENDC}")
    wave_merged, flux_merged = read_baches_merged(file_with_merged_spectrum)

    ###
    #   Normalize flux
    #
    if normalize:
        print(
            f"{Bcolors.BOLD}   Normalize merged MIDAS spectrum{Bcolors.ENDC}"
        )
        if norm_version == 'specutils_continuum':
            wave_merged, flux_merged = normalize_spectrum_interval(
                wave_merged,
                flux_merged,
                median_window,
                porder,
            )
            wave_merged = wave_merged.value
            flux_merged = flux_merged.value
            continuum = None
        elif norm_version == 'median_max_window':
            wave_merged = wave_merged.value
            flux_merged = flux_merged.value
            flux_merged, continuum = normalize_spectrum_fabian(
                wave_merged,
                flux_merged,
                object_name=object_name,
            )
        else:
            print(
                f"{Bcolors.FAIL}   Normalize method not known. Check variable:"
                f" 'norm_version'. Skipping normalization!{Bcolors.ENDC}"
            )
            continuum = None
    else:
        wave_merged = wave_merged.value
        flux_merged = flux_merged.value
        continuum = None

    ###
    #   Correct for radial velocity -> 1) calculate barycentric velocity
    #                                     correction
    #                                  2) add radial velocity
    #
    velocity_correction = radial_velocity
    if correct_bary:
        try:
            bvc = bary_correction(file_with_merged_spectrum)
        except RuntimeError:
            print(
                f"{Bcolors.WARNING}   Barycentric velocity correction could "
                f"not be determined. Assume 0 km/s.{Bcolors.ENDC}"
            )
            bvc = 0.
        velocity_correction = radial_velocity - bvc

    wave_merged = correct_for_radial_velocity(
        wave_merged,
        flux_merged,
        velocity_correction,
        object_name,
        # plot=True,
        merge_type='MIDAS',
    )

    ###
    #   Plot merged data in one plot
    #
    print(f"{Bcolors.BOLD}   Plot spectra{Bcolors.ENDC}")
    plot_merged(
        wave_merged,
        flux_merged,
        object_name,
        normalize,
        merge_type='MIDAS',
    )

    ###
    #   Plot merged data in individual panels
    #
    if panel_version == 'old':
        plot_panels(
            wave_merged,
            flux_merged,
            object_name,
            normalize,
            line_file,
            panel_wave_range,
            merge_type='MIDAS',
        )
    elif panel_version == 'default':
        #   Generate ident lines
        lines = generate_lines_from_file(
            wave_merged,
            flux_merged,
            ions,
            line_file=line_file,
            percentage_line_flux_must_be_below_continuum=percentage_line_flux_must_be_below_continuum,
        )

        plot_panels_fabian(
            wave_merged,
            flux_merged,
            object_name,
            normalize,
            {**lines, **manual_lines},
            n_panels=n_panels,
            merge_type='MIDAS',
            radial_velocity=radial_velocity,
            continuum=continuum,
            n_panels_per_page=n_panels_per_page,
        )
    else:
        print(
            f"{Bcolors.FAIL}   Version for panel version not known. Check "
            f"variable: 'panel_version'. Skipping panel plot!{Bcolors.ENDC}"
        )

    ########################################################################
    #                           Individual orders                          #
    ########################################################################
    if individual_orders:
        print(f"{Bcolors.BOLD}Processing individual orders{Bcolors.ENDC}")

        ###
        #   Load spectra
        #
        #   Select FIT or fit
        file_with_orders = check_file_name(file_with_orders)

        #   Open fits file and read the data and header section
        print(f"{Bcolors.BOLD}   Read files ...{Bcolors.ENDC}.")
        wave_list, flux_list = read_baches(file_with_orders)

        if normalize and norm_individual:
            ###
            #   Normalize individual orders and merge afterwards
            #
            print(
                f"{Bcolors.BOLD}   Normalize individual "
                f"orders...{Bcolors.ENDC}"
            )

            orders = []
            for i, wave in enumerate(wave_list):
                orders.append(
                    Spectrum1D(spectral_axis=wave, flux=flux_list[i])
                )

            #   Normalize & merge spectra
            merged_wave, merged_flux = normalize_and_merge_spectra(
                orders,
                median_window=median_window,
                normalization_order=porder,
            )
        else:
            ###
            #   Merge orders
            #
            if debug_plot:
                #   Prepare plot
                fig = plt.figure(figsize=(18, 9))
                plt.xlabel(r'Wavelength $\lambda\,[\AA]$')
                plt.ylabel('Relative flux')
                plt.grid(visible=True, axis='y')

            #   Process each order
            merged_wave, merged_flux = merge_spectra(
                wave_list,
                flux_list,
                trim_value_start_order=trim_value_start_order,
                debug_plot=debug_plot,
            )

            if debug_plot:
                #   Plot merged spectrum
                plt.step(merged_wave, merged_flux, ls='--')

                plt.show()
                plt.close()

        ###
        #   Normalize merged spectrum
        #
        if normalize and not norm_individual:
            print(
                f"{Bcolors.BOLD}   Normalize merged spectrum{Bcolors.ENDC}"
            )
            if norm_version == 'specutils_continuum':
                merged_wave, merged_flux = normalize_spectrum_interval(
                    merged_wave,
                    merged_flux,
                    median_window,
                    porder,
                )
                merged_wave = merged_wave.value
                merged_flux = merged_flux.value
                continuum = None
            elif norm_version == 'median_max_window':
                merged_wave = merged_wave.value
                merged_flux = merged_flux.value
                merged_flux, continuum = normalize_spectrum_fabian(
                    merged_wave,
                    merged_flux,
                    object_name=object_name,
                )
            else:
                print(
                    f"{Bcolors.FAIL}   Normalize method not known. "
                    f"Check variable: 'norm_version'. Skipping "
                    f"normalization!{Bcolors.ENDC}"
                )
        else:
            merged_wave = merged_wave.value
            merged_flux = merged_flux.value
            continuum = None

        ###
        #   Correct for radial velocity
        #
        merged_wave = correct_for_radial_velocity(
            merged_wave,
            merged_flux,
            velocity_correction,
            object_name,
            # plot=True,
            merge_type='PYTHON',
        )

        ###
        #   Plot merged data in one plot
        #
        print(f"{Bcolors.BOLD}   Plot spectra{Bcolors.ENDC}")
        plot_merged(
            merged_wave,
            merged_flux,
            object_name,
            normalize,
            merge_type='PYTHON',
        )

        ###
        #   Plot merged data in individual panels
        #
        if panel_version == 'old':
            plot_panels(
                merged_wave,
                merged_flux,
                object_name,
                normalize,
                line_file,
                panel_wave_range,
                merge_type='PYTHON',
            )
        elif panel_version == 'default':
            #   Generate ident lines
            lines = generate_lines_from_file(
                merged_wave,
                merged_flux,
                ions,
                line_file=line_file,
                percentage_line_flux_must_be_below_continuum=percentage_line_flux_must_be_below_continuum,
            )

            plot_panels_fabian(
                merged_wave,
                merged_flux,
                object_name,
                normalize,
                {**lines, **manual_lines},
                n_panels=n_panels,
                merge_type='PYTHON',
                radial_velocity=radial_velocity,
                continuum=continuum,
                n_panels_per_page=n_panels_per_page,
            )
        else:
            print(
                f"{Bcolors.FAIL}   Version for panel version not known. "
                f"Check variable: 'panel_version'. Skipping panel "
                f"plot!{Bcolors.ENDC}"
            )

    print(f'{Bcolors.OKGREEN}DONE{Bcolors.ENDC}')
