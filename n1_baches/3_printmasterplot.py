#! /usr/bin/python3
# -*- coding: utf-8 -*-

############################################################################
#             Configuration: modify the file in this section               #
############################################################################

#   Name of file with individual orders
File_orders = "master_spectrum_wr.fit"

#   Name of file with merged spectrum
File_merged = "master_spectrum_wrm.fit"

#   Name of the object
objectname = "star"

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
line_lower_continuum = 3.

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
n_panel = 5

#   Number of panels into which the spectrum should be split
#   (panel_version=default)
divinto = 21

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
#   be normalized afterwards.
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
trim_value = 0
# trim_value = 400


###
#   Line identifications
#
#   File containing line identifications
# lineFile = ""
# lineFile = "absorption_lines.dat"
lineFile = "/home/pollux/reduction_scripts_students/n1_baches/atomic_lines.tsv"

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


def consecutive(data, stepsize=1):
    """
        Find consecutive elements in a numpy array

        Idee       : stackoverflow.com/questions/7352684/how-to-find-the-groups-of-consecutive-elements-in-a-numpy-array
        ----

        Parameters
        ----------
        data       : `numpy.ndarray`
            Data array

        stepsize   : `integer`, optional
            Step size between the consecutive elements
    """
    return np.split(data, np.where(np.diff(data) != stepsize)[0] + 1)


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
        obs_time = Time(date, format='jd', scale='utc')
    except RuntimeError:
        date = header['DATE-OBS']
        obs_time = Time(date, format='fits', scale='utc')

    ra = header['OBJCTRA']
    dec = header['OBJCTDEC']
    observatory_latitude = header.get('SITELAT', 52.4)
    observatory_longitude = header.get('SITELONG', 12.96)
    observatory_height_above_MSL = header.get('SITEHIGH', 40.)

    #   Define Location
    obs_site = EarthLocation(
        lat=observatory_latitude,
        lon=observatory_longitude,
        height=observatory_height_above_MSL,
    )

    #   Calculate barycentric velocity correction with the help of a
    #   SkyCoord object
    obs_coord = SkyCoord(
        ra=ra,
        dec=dec,
        unit=(u.hourangle, u.deg),
        frame='icrs',
        obstime=obs_time,
        location=obs_site
    )

    bvc = obs_coord.radial_velocity_correction(kind='barycentric').to(u.km / u.s)

    return bvc.value


def correct_for_rv(wl, flx, rv, name_obj, plot=False, mtype=None):
    """
        Correct wave length array for radial velocity

        Parameters
        ----------
        wl              : `numpy.ndarray`
            Wavelength data

        flx             : `numpy.ndarray`
            Flux data

        rv              : `float`
            Radial velocity

        name_obj        : `string`
            Name of the object

        plot            : `boolean`, optional
            If `True` the plot will be shown
            Default is ``False``.

        mtype           : `string`, optional
            String that characterize the order merging procedure. It that will
            be part of the file names.
            Default is ``None``.

        Returns
        -------
        corr_wl         : `numpy.ndarray`
            Redial velocity corrected wavelength
    """
    #   Correct wave length data
    if rv:
        corr_wl = wl / (np.sqrt((c + rv * 1000) / (c - rv * 1000)))
    else:
        corr_wl = wl

    #   Get plot limits -> use line wavelength + a margin of 20AA
    id_x_min = np.argmin(np.abs(corr_wl - (5889.95 - 20)))
    id_x_max = np.argmin(np.abs(corr_wl - (5895.92 + 20)))

    plt.plot(
        corr_wl[id_x_min:id_x_max],
        flx[id_x_min:id_x_max],
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
    if mtype is not None:
        plt.savefig(f"output/Vicinity_sodium_{name_obj}_{mtype}.pdf")
    else:
        plt.savefig(f"output/Vicinity_sodium_{name_obj}.pdf")
    plt.close()

    return corr_wl


def normfunc(p, f1, f2):
    """
        Function to minimize flux difference between spectra

        Goal        : Calculate the number of wavelength points for which the
        ----          flux difference is below a certain threshold, so this can
                      be maximized.

        Parameters
        ----------
        p           : `float`
            Initial value for the flux scaling between spectra

        f1          : `numpy.ndarray`
            Flux of spectrum 1

        f2          : `numpy.ndarray`
            Flux of spectrum 2
    """
    #   Set flux threshold
    limit = np.abs(0.02 * np.nanmedian(f1))

    #   Calculate points where flux differences is below the threshold
    good_points = (np.abs(f1 - f2 * p) < limit).sum()

    #   Sanitize number of points, since we will return the inverse
    if not good_points:
        good_points = 0.00001

    return 1. / good_points


def norm_spectrum_interval(wave, flux, median_window, porder):
    """
        Normalize a spectrum. If the normalization is not sufficient, the
        spectrum will be split in segments and the normalization will
        be repeated for those segments. The segments will be merged afterwards.

        Parameters
        ----------
        wave            : `numpy.ndarray`
            Array with wavelength data

        flux            : `numpy.ndarray`
            Array with flux data

        median_window   : `int`
            Window in Pixel used in median smoothing

        porder          : `integer`
            Order of the polynomial that is used for normalization

        Returns
        -------
        wave            : `numpy.ndarray`
            Array with normalized wavelength data

        flux            : `numpy.ndarray`
            Array with normalized  flux data
    """
    #   Create Spectrum1D objects
    spec = Spectrum1D(spectral_axis=wave, flux=flux)

    #   Normalize the merged spectrum
    spec, std = norm_spectrum(spec, order=porder)

    #   Split spectrum in segments, if standard deviation is too high
    if std > 0.05:
        nsegment = 10
        nwave = len(wave)
        step = int(nwave / nsegment)
        segments = []

        #   Loop over segments
        i_old = 0
        for i in range(step, step * nsegment, step):
            #   Cut segments and add overlay range to the
            #   segments, so that the normalization afterburner
            #   can take effect
            overlap = int(step * 0.15)
            if i == step:
                flux_seg = flux[i_old:i + overlap]
                wave_seg = wave[i_old:i + overlap]
            elif i == nsegment - 1:
                flux_seg = flux[i_old - overlap:]
                wave_seg = wave[i_old - overlap:]
            else:
                flux_seg = flux[i_old - overlap:i + overlap]
                wave_seg = wave[i_old - overlap:i + overlap]
            i_old = i

            #   Create Spectrum1D objects for the segments
            segments.append(Spectrum1D(spectral_axis=wave_seg, flux=flux_seg))

        #   Normalize & merge spectra
        wave, flux = norm_merge_spectra(
            segments,
            median_window=median_window,
            order=porder,
        )

    return wave, flux


def norm_spectrum(spec, median_window=61, order=9):
    """
        Normalize a spectrum

        Parameters
        ----------
        spec            : `specutils.Spectrum1D`
            Spectrum to normalize

        median_window   : `int`, optional
            Window in Pixel used in median smoothing
            Default is ``61``.

        order           : `int`, optional
            Order of the polynomial used to find the continuum
            Default is ``9``.

        Returns
        -------
        norm_spec       : `specutils.Spectrum1D`
            Normalized spectrum
    """
    #   Regions that should not be used for continuum estimation,
    #   such as broad atmospheric absorption bands
    exclude_regions = [
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
    _cont = fit_generic_continuum(
        spec,
        model=models.Chebyshev1D(order),
        fitter=fitting.LinearLSQFitter(),
        median_window=median_window,
        exclude_regions=exclude_regions,
    )(spec.spectral_axis)

    #   Normalize spectrum
    norm_spec = spec / _cont

    #   Sigma clip the normalized spectrum to rm spectral lines
    clip_flux = sigma_clip(
        norm_spec.flux.value,
        sigma_lower=1.25,
        sigma_upper=3.,
        # sigma_upper=4.,
        axis=0,
        grow=1.,
    )

    #   Calculate mask
    mask = np.invert(clip_flux.recordmask)

    #   Make new spectrum
    spec_mask = Spectrum1D(
        spectral_axis=spec.spectral_axis[mask],
        flux=spec.flux[mask],
    )

    # Determine new continuum
    _cont = fit_generic_continuum(
        spec_mask,
        model=models.Chebyshev1D(order),
        fitter=fitting.LinearLSQFitter(),
        median_window=median_window,
        exclude_regions=exclude_regions,
    )(spec.spectral_axis)

    #   Normalize spectrum again
    norm_spec = spec / _cont

    return norm_spec, mad_std(norm_spec.flux)


def norm_merge_spectra(spectra, median_window=61, order=9):
    """
        Normalize spectra and merge them afterwards

        Parameters
        ----------
        spectra         : `list` of `specutils.Spectrum1D`
            Spectra to normalize and merge

        median_window   : `int`, optional
            Window in Pixel used in median smoothing
            Default is ``61``.

        order           : `int`, optional
            Polynomial order used for normalization
            Default is ``9``.

        Returns
        -------
                        :
            Normalized and merged spectrum
    """
    #   Normalize spectra
    norm_spec = []
    for spec in spectra:
        norm_spec.append(
            norm_spectrum(
                spec,
                median_window=median_window,
                order=order
            )[0]
        )

    #   Merge spectra
    return merge_norm(norm_spec)


def merge_norm(spec_list):
    """
        Merge normalized spectra

        Parameters
        ----------
        spec_list       : `list` of `specutils.Spectrum1D`
            List with normalized spectra

        Returns
        -------
        mwave           : `list` of `float`
            Common wave length scale

        mflux          : `list` of `float`
            Merged flux data

    """
    #   Extract wavelength and flux data
    wave = []
    flux = []
    for spec in spec_list:
        wave.append(spec.spectral_axis)
        flux.append(spec.flux)

    #   Build common wave length scale
    mwave = np.sort(np.concatenate(wave))

    #   Prepare list for flux
    flux_new = []

    #   Get flux unit
    flux_u = flux[0].unit

    #   Loop over all spectra
    for i, w in enumerate(wave):
        #   Remove wavelength duplicates from the input arrays
        _u, indices = np.unique(w, return_index=True)
        w = w[indices]
        flux[i] = flux[i][indices]

        #   Interpolate flux on new wavelength grid
        f = interpolate.interp1d(
            w,
            flux[i],
            kind='cubic',
            bounds_error=False,
        )
        flux_new.append(f(mwave))

    #   Normalization afterburner
    #   -> attempts to improve normalization in the overlap regions
    flux_new = after_norm(flux_new)

    #   Merge flux
    mflux = np.nanmedian(flux_new, axis=0)

    return mwave, mflux * flux_u


def after_norm(flux_list):
    """
        Afterburner that attempts to correct spectra (especially echelle
        orders) that are not well normalized before merging.
            -> only works if there in an overlap region
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
        f_diff = flux_list[i - 1] - flux_list[i]

        #   Check if overlap region exists
        #   -> if not do nothing
        if ~np.all(np.isnan(f_diff)):
            #   Determine where flux is negative
            lower = np.argwhere(f_diff < -0.01)

            #   Find consecutive elements -> build groups
            group_lower = consecutive(lower.flatten())

            #   Loop over groups
            for low in group_lower:
                #   Restrict to groups with more than 5 elements
                if len(low) > 5:
                    #   Check if flux of the first spectrum
                    #   is negative in this area
                    #   -> if yes replace
                    #   -> otherwise replace flux in the second spectrum
                    if np.nanmedian(flux_list[i - 1][low] - 1) < -0.05:
                        flux_list[i - 1][low] = flux_list[i][low]
                    else:
                        flux_list[i][low] = flux_list[i - 1][low]

            #   Determine where flux is positive
            higher = np.argwhere(f_diff > 0.01)

            #   Find consecutive elements -> build groups
            group_higher = consecutive(higher.flatten())

            #   Loop over groups
            for high in group_higher:
                #   Restrict to groups with more than 5 elements
                if len(high) > 5:
                    #   Check if flux of the second spectrum
                    #   is negative in this area
                    #   -> if yes replace
                    #   -> otherwise replace flux in the first spectrum
                    if np.nanmedian(flux_list[i][high] - 1) < -0.05:
                        flux_list[i][high] = flux_list[i - 1][high]
                    else:
                        flux_list[i - 1][high] = flux_list[i][high]

    return flux_list


def merge_spectra_intersec(flux_1, flux_2, debug_plot=False, wave=None,
                           i=None):
    """
        Combines two spectra whose flux intersects

        Goal        : Tries to avoid or at least reduce jumps, while merging
        ----

        Parameters
        ----------
        flux_1      : `numpy.ndarray`
            Flux of the first spectrum

        flux_2      : `numpy.ndarray`
            Flux of the second spectrum

        debug_plot  : `boolean`, optional
            If ``True`` a the intersection range will be marked on the
            debug plot
            Default is ``False``.

        wave        : `numpy.ndarray`, optional
            Wavelength array
            Default is ``None``.

        i           : `integer`, optional
            ID of the current "intersection"
            Default is ``None``.

        Returns
        -------
        flux_1      : `numpy.ndarray`
            Modified flux of the first spectrum

        flux_2      : `numpy.ndarray`
            Modified flux of the second spectrum
    """
    #   Calculate flux difference
    f_diff = flux_1 - flux_2

    #   Check if overlap exists?
    if ~np.all(np.isnan(f_diff)):
        #   Calculate range where the flux difference < 10% of maximum of
        #   the flux difference to identify the range where flux intersects
        id_f_x = np.argwhere(
            np.absolute(f_diff) / np.nanmax(np.absolute(f_diff)) <= 0.1
        )

        #   Return is no intersection was identified
        if len(id_f_x) == 0:
            return flux_1, flux_2

        #   Find center of id_f_x
        len_f_x = len(id_f_x)
        id_f_x = id_f_x[int(len_f_x / 2)][0]

        #   Median flux around 10% of id_f_x
        #   (plus 1 to ensure that range is not 0)
        x_len = int(len_f_x * 0.05) + 1
        flux_x = np.median(flux_1[id_f_x - x_len:id_f_x + x_len])
        if flux_x == 0.:
            flux_x = 1.

        #   Cut flux range to the overlap range
        f_diff_cut = f_diff[~np.isnan(f_diff)]

        #   First and last flux point
        f_diff_s = f_diff_cut[0]
        f_diff_e = f_diff_cut[-1]

        #   Find index of the above
        id_f_s = np.nanargmin(np.absolute(f_diff - f_diff_s))
        id_f_e = np.nanargmin(np.absolute(f_diff - f_diff_e))

        #   Plot markers for overlapping region
        if debug_plot:
            plt.plot(
                [wave[id_f_s].value, wave[id_f_s].value],
                [0.9 * flux_1[id_f_s], 1.1 * flux_1[id_f_s]],
                color='b',
                linestyle='--',
            )
            plt.plot(
                [wave[id_f_e].value, wave[id_f_e].value],
                [0.9 * flux_1[id_f_e], 1.1 * flux_1[id_f_e]],
                color='b',
                linestyle='--',
            )
            plt.plot(
                [wave[id_f_x].value, wave[id_f_x].value],
                [0.9 * flux_1[id_f_x], 1.1 * flux_1[id_f_x]],
                color='r',
                linestyle='--',
            )
            plt.text(
                wave[id_f_x].value,
                1.2 * flux_1[id_f_x],
                ' ' + str(i) + ' ',
                # rotation=90,
                ha='center',
                va='top',
                fontdict=font,
            )

        #   Calculate 3% of the length of the overlap range
        three_diff = int(len(f_diff_cut) * 0.03)

        #   Ensure that the above is at least 1
        if three_diff == 0:
            three_diff = 1

        #   Median flux of the first and last 3% in terms of bins
        f_diff_s_med = np.median(f_diff_cut[0:three_diff])
        f_diff_e_med = np.median(f_diff_cut[three_diff * -1:])

        #   Check if flux difference stars negative and ends positive
        #   and is grater than 3% of the median flux
        #   -> if yes, use flux of the other in the respective area
        # if (f_diff_s_med / flux_x < -0.03 and f_diff_e_med / flux_x > 0.03 or
        # f_diff_s_med / flux_x > 0.03 and f_diff_e_med / flux_x < -0.03):
        # flux_2[id_f_s:id_f_x] = flux_1[id_f_s:id_f_x]
        # flux_1[id_f_x:id_f_e] = flux_2[id_f_x:id_f_e]

        #   If the flux difference is larger than 3% at one edge, remove
        #   the overlapping flux ranges, since the following merging process
        #   would in this case lead to jumps in the merged spectrum
        # if (np.abs(f_diff_s_med / flux_x) < 0.3 or np.abs(f_diff_e_med / flux_x) < 0.3):
        if (np.abs(f_diff_s_med / flux_x) > 0.03 or
                np.abs(f_diff_e_med / flux_x) > 0.03):
            flux_2[id_f_s:id_f_x] = np.nan
            flux_1[id_f_x:id_f_e] = np.nan

    return flux_1, flux_2


def norm_two_spectra(flux_1, flux_2, p0):
    """
        Normalize and adjust flux of two spectra

        Parameters
        ----------
        flux_1          : `numpy.ndarray`
            Flux of first spectrum

        flux_2          : `numpy.ndarray`
            Flux of second spectrum

        p0              : `float`
            Initial value for the minimization algorithm used to estimated
            the flux difference between `flux_1` and `flux_2`. Usually the
            normalization factor of previous iteration/order

        Returns
        -------

                        : `numpy.ndarray`
            Flux of second spectrum (`flux_2`) scaled to the first
            one (`flux_1`)

                        : `float`
            Scaling factor of `flux_2`

        Idea            : https://stackoverflow.com/questions/13846213/create-composite-spectrum-from-two-unnormalized-spectra
        ----
    """

    #   Calculate normalization factor between `flux_1` and `flux_2`, using
    #   scipy.optimize and a normalization function `normfunc`
    min_fit = optimize.basinhopping(
        # min_fit = optimize.minimize_scalar(
        normfunc,
        p0,
        1000,
        # args=(flux_1, flux_2),
        minimizer_kwargs={'args': (flux_1, flux_2)}
        # )
    )
    #   If the minimization algorithm fails, use the normalization factor of
    #   previous order
    if min_fit.success:
        norm_factor = min_fit.x
    else:
        norm_factor = p0

    #   Adjust flux and return
    return flux_2 * norm_factor, norm_factor


def merge_spectra(wave, flux, trim_value=400, debug_plot=False):
    """
        Resample spectra on a common wavelength scale and merge those

        Parameters
        ----------
        wave            : `list` of `numpy.ndarray`'s
            List of the wavelength ranges

        flux            : `list` of `numpy.ndarray`'s
            List of the flux ranges corresponding to the wavelength ranges
            in `wave` that should be merged

        trim_value      : `integer`, optional
            The number of wavelength bins that should be removed from the start
            of the orders.
            Default is ``400``.

        debug_plot  : `boolean`, optional
            If ``True`` a the intersection range will be marked on the
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
    new_wave = np.unique(np.sort(np.concatenate(wave)))

    #   Prepare list for flux
    flux_new = []

    #   Initial normalization value
    p0 = 1.

    #   Get flux unit
    flux_u = flux[0].unit

    print(f"{Bcolors.BOLD}   Match order fluxes ...{Bcolors.ENDC}")
    for i, w in enumerate(wave):
        #   Calculate positions where flux is 0. and rm those
        ind = np.argwhere(flux[i] == 0.)

        wave_clean = np.delete(w, ind)
        flux_clean = np.delete(flux[i], ind)

        #   Trim start of the spectral ranges to improve normalization
        #   This proves to be useful in some cases and does not do harm in
        #   the other.
        if i == 0:
            wave_clean = wave_clean[600:]
            flux_clean = flux_clean[600:]
        if i > 0:
            wave_clean = wave_clean[trim_value:]
            flux_clean = flux_clean[trim_value:]

        #   Interpolate flux on common wavelength scale
        f = interpolate.interp1d(
            wave_clean,
            flux_clean,
            kind='cubic',
            bounds_error=False,
        )
        flux_new.append(f(new_wave))

        if i > 0:
            print(f"   Adjust flux of order {i} to {i + 1}\r", end="")
            #   Adjust flux of the individual spectra
            flux_new[i], p0 = norm_two_spectra(
                flux_new[i - 1],
                flux_new[i],
                p0,
            )
        if debug_plot:
            plt.step(new_wave, flux_new[i])

    #   Improve overlapping edges before merging
    for i, w in enumerate(wave):
        if i > 0:
            print(
                f"   Improve overlapping edges for order {i} and {i + 1}\r",
                end=""
            )
            #   Manipulate spectra if they "intersect"
            flux_new[i - 1], flux_new[i] = merge_spectra_intersec(
                flux_new[i - 1],
                flux_new[i],
                debug_plot=debug_plot,
                wave=new_wave,
                i=i,
            )
    print("                                                        \r", end="")

    #   Merge flux and remove residuals NANS
    print(f"{Bcolors.BOLD}   Merge orders... {Bcolors.ENDC}")
    mflux = np.nanmedian(flux_new, axis=0)
    ind = np.argwhere(np.isnan(mflux))

    return np.delete(new_wave, ind), np.delete(mflux, ind) * flux_u


def normalize_spectrum_fabian(wl, flx, name="", cutlines=[], cutlim=[]):
    """
        Normalize a spectrum to its continuum

        Parameters
        ----------
        wl              : `numpy.ndarray`
            Wavelength data

        flx             : `numpy.ndarray`
            Flux data

        name            : `string`, optional
            Name of the object
            Default is an empty string.

        cutlines        : `list`, optional
            Line wavelength that should be skipped during continuum
            determination
            Default is an empty list.

        cutlim          : `list`, optional
            Wavelength limits that should be applied around each line that
            should be skipped
            Default is an empty list.

        Returns
        -------
        n_flx           : `numpy.ndarray`
            Normalized flux data

        cont            : `numpy.ndarray`
            Continuum data points
    """
    #   Calculate median of the step size between wavelength points
    wl_step = np.median(np.diff(wl))

    med_win_size = 0.5
    max_win_size = 12
    #   Calculate number of elements in the maximum and median window
    true_med_size = math.floor(med_win_size / wl_step)
    true_max_size = math.floor(max_win_size / wl_step)

    if true_max_size == 0 or true_med_size == 0:
        raise ValueError(
            "Medium/Maximum window sizes need to be bigger than a "
            "wavelength step!"
        )

    #   Apply median and maximum filter with the corresponding window sizes
    flx_for_interpol = median_filter(flx, size=true_med_size)
    flx_for_interpol = maximum_filter(flx_for_interpol, size=true_max_size)

    #   Cut spectral lines that would spoil the continuum fit such as broad
    #   hydrogen lines
    cutmask = np.ones(np.shape(wl), dtype=bool)
    for i, l in enumerate(cutlines):
        line_mask = np.logical_or(wl > l + cutlim[i], wl < l - cutlim[i])
        cutmask = np.logical_and(cutmask, line_mask)

    flx_for_interpol = flx_for_interpol[cutmask]
    wl_for_interpol = wl[cutmask]

    #   Interpolate to find the continuum
    norm_fit = InterpolatedUnivariateSpline(
        wl_for_interpol,
        flx_for_interpol,
        k=2,
    )

    #   Normalize spectrum
    cont = norm_fit(wl)
    n_flx = flx / cont

    if name != "":
        np.savetxt(
            f"flux_normalized_{name}.csv",
            np.transpose((wl, n_flx)),
            delimiter=",",
        )

    return n_flx, cont


def read_baches(File_orders):
    """
        Read baches FITS files with individual orders created by MIDAS

        Parameters
        ----------
        File_orders        : `string`
            Name of the file to read

        Returns
        -------
        wave_list       : `list of `numpy.ndarray`
            List with the wavelength data of the individual orders

        flux_list       : `list of `numpy.ndarray`
            List with the flux data of the individual orders
    """
    File = fits.open(File_orders)
    FileData = File[0].data

    FileHeader = File[0].header
    History = FileHeader['HISTORY']
    # ref_pixel        = FileHeader['CRPIX1']
    # coord_ref_pixel  = FileHeader['CRVAL1']
    wave_per_pixel = FileHeader['CDELT1']

    #   Extract startpoints from HISTORY section of the FITS Header
    k = 1
    startpoints = []
    for line in History:
        liste = line.split()
        if not liste:
            k = 1
            continue
        if k == 0:
            startpoints = startpoints + liste
        if k == 1:
            if 'WSTART' not in liste[0]:
                continue
            else:
                k = 0
    startpoints = np.array(startpoints, dtype=float)

    #   Convert fluxes to numpy array and calculate wavelength points for
    #   all flux values
    fluxes = np.array(FileData, dtype=float)
    waves = np.ones(fluxes.shape) * np.arange(0, fluxes.shape[1]) * wave_per_pixel
    waves = waves + startpoints.reshape(startpoints.size, 1)

    #   Limit to flux values != 0
    wave_list = []
    flux_list = []
    for i, flux in enumerate(fluxes):
        mask_id = np.nonzero(flux)
        wave_list.append(waves[i, mask_id][0] * u.AA)
        flux_list.append(flux[mask_id] * u.adu)

    return wave_list, flux_list


def read_baches_merged(File_orders):
    """
        Read baches FITS file with spectrum merged by MIDAS

        Parameters
        ----------
        File_orders        : `string`
            Name of the file to read

        Returns
        -------
        wave_merged     : `numpy.ndarray`
            Wavelength data of the spectrum merged by MIDAS

        flux_merged     : `numpy.ndarray`
            Flux data of the spectrum merged by MIDAS
    """
    File_merged = fits.open(File_orders)
    FileData_merged = File_merged[0].data

    FileHeader_merged = File_merged[0].header
    # ref_pixel_merged       = FileHeader_merged['CRPIX1']
    coord_ref_pixel_merged = FileHeader_merged['CRVAL1']
    wave_per_pixel_merged = FileHeader_merged['CDELT1']

    #   Number of lines in the file
    FileDataLines_merged = len(FileData_merged)
    N2 = FileDataLines_merged

    wave_merged = []
    for i in range(0, FileDataLines_merged):
        wave_merged.append(coord_ref_pixel_merged + i * wave_per_pixel_merged)

    scut = 200
    ecut = 200
    wave_merged = wave_merged[scut:-ecut]
    flux_merged = FileData_merged[scut:-ecut]

    return wave_merged * u.AA, flux_merged * u.adu


def check_file_name(File):
    """
        Check that the file has the correct ending

        Parameters
        ----------
        File            : `string`
            Path to the file

        Returns
        -------
        File            : `string`
            Modified path to the file
    """
    if os.path.isfile(File) is False:
        if File.find(".fit") > 0:
            File = File.replace(".fit", ".FIT")
            print(f"change .fit to .FIT in {File}")
        elif ".FIT" in File:
            File = File.replace(".FIT", ".fit")
            print(f"change .FIT to .fit in {File}")
        else:
            print(
                f"{Bcolors.FAIL}!!!!! Check the right spelling of the file "
                f"extension and the filepath !!!!!{Bcolors.ENDC}"
            )

    return File


def add_idents(wave, flux, lineFile, minflux=None, maxflux=None):
    """
        Add idents to the plot

        Parameters
        ----------
        wave            : `numpy.ndarray`
            Wavelength array

        flux            : `numpy.ndarray`
            Flux array

        minflux         : `float`, optional
            Minimum plot range on the Y axis
            Default is ``None``.

        maxflux         : `float`, optional
            Maximum plot range on the Y axis
            Default is ``None``.

        lineFile        : `string`
            Name of the ident file
    """
    #   Defining plot range
    if minflux is None:
        minflux = min(flux)
        minflux = max(0, minflux)
    if maxflux is None:
        maxflux = max(flux)

    yoffset = (maxflux - minflux) * 0.01
    xoffset = (max(wave) - min(wave)) * 0.01

    plotminimum = minflux - yoffset
    plotmaximum = maxflux + yoffset

    #   Setting plotpositions for ident lines
    plotheigth = plotmaximum - plotminimum
    plotmiddleplot = (plotminimum + plotmaximum) / 2.0
    plotupperplot = plotminimum + 0.90 * plotheigth
    plotlowerplot = plotminimum + 0.10 * plotheigth
    plotuppercut1 = plotminimum + 0.80 * plotheigth
    plotlowercut1 = plotminimum + 0.20 * plotheigth
    plotuppercut2 = plotminimum + 0.78 * plotheigth
    plotlowercut2 = plotminimum + 0.22 * plotheigth

    #   Interpolate on data to find point for ident
    f2 = interpolate.interp1d(wave, flux)

    #   Open ident file
    try:
        lines = open(lineFile, "r")
    except RuntimeError as e:
        print(
            f"{Bcolors.FAIL}   Line file not found. Check variable "
            f"'lineFile'. Specified was {lineFile}. {Bcolors.ENDC}"
        )
        raise e

    #   Plot idents to figure
    for line in lines:
        liste = line.split()

        if len(liste) == 1:
            print(
                f"{Bcolors.WARNING}   [WARNING] Broken identification found "
                f"as '{line}', must consist of [wavelength(s) + name]. I will "
                f"skip this one.{Bcolors.ENDC}"
            )
            continue
        try:
            float(liste[0])
        except ValueError:
            print(
                f"{Bcolors.WARNING}   [WARNING] Broken identification found "
                f"as '{line}', first entry not a number. I will skip this "
                f"one.{Bcolors.ENDC}"
            )
            continue

        #   Single ident
        if len(liste) == 2:
            ident_line = liste

            #   Only plot if in range
            if float(ident_line[0]) >= min(wave) and float(ident_line[0]) <= max(wave):
                #   Plot ident point in upper or lower half of figure
                if f2(ident_line[0]) >= plotmiddleplot:
                    plt.plot(
                        [float(ident_line[0]), float(ident_line[0])],
                        [f2(float(ident_line[0])), plotlowercut2],
                        color='r',
                        linestyle='-',
                    )
                    plt.text(
                        float(ident_line[0]),
                        plotlowercut2,
                        ' ' + ident_line[1] + ' ',
                        rotation=90,
                        ha='center',
                        va='top',
                        fontdict=font,
                    )
                else:
                    plt.plot(
                        [float(ident_line[0]), float(ident_line[0])],
                        [plotuppercut2, f2(float(ident_line[0]))],
                        color='r',
                        linestyle='-',
                    )
                    plt.text(
                        float(ident_line[0]),
                        plotuppercut2,
                        ' ' + ident_line[1] + ' ',
                        rotation=90,
                        ha='center',
                        va='bottom',
                        fontdict=font,
                    )

        #   Multi ident
        if len(liste) > 2:
            points = []
            pointsFlux = []
            for i in liste[:-1]:
                points.append(float(i))
            pointcenter = sum(points) / float(len(points))
            ident_name = str(liste[-1:])
            if max(points) <= max(wave) and min(points) >= min(wave):
                for i in liste[:-1]:
                    pointsFlux.append(f2(float(i)))
                pointFluxCenter = sum(pointsFlux) / float(len(points))

                if pointFluxCenter <= plotmiddleplot:
                    plt.plot(
                        [pointcenter, pointcenter],
                        [plotupperplot, plotuppercut1],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        pointcenter,
                        plotupperplot,
                        ' ' + ident_name[2:-2] + ' ',
                        rotation=90,
                        ha='center',
                        va='bottom',
                        fontdict=font,
                    )
                    for element in points:
                        plt.plot(
                            [element, element],
                            [f2(element), plotuppercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )
                        plt.plot(
                            [pointcenter, element],
                            [plotuppercut1, plotuppercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )

                if pointFluxCenter > plotmiddleplot:
                    plt.plot(
                        [pointcenter, pointcenter],
                        [plotlowerplot, plotlowercut1],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        pointcenter,
                        plotlowerplot,
                        ' ' + ident_name[2:-2] + ' ',
                        rotation=90,
                        ha='center',
                        va='top',
                        fontdict=font,
                    )
                    for element in points:
                        plt.plot(
                            [element, element],
                            [f2(element), plotlowercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )
                        plt.plot(
                            [pointcenter, element],
                            [plotlowercut1, plotlowercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.0,
                        )
    lines.close()


def add_idents_fabian(lines, xlim, axs):
    """
        Add idents to the plot

        Parameters
        ----------
        lines                   : `dictionary`
            Line information: key=line identifier, value[0]=wavelength,
                                                   value[1]=alignment parameter

        xlim                    : `tupel`
            Limits for the plot in X direction

        axs                     : `matplotlib.pyplot.subplots`
            Plot to which the idents should be added.
    """
    for linestr, linelocs in lines.items():
        alllinelocs = [l[0] for l in linelocs]
        #   Restrict to line within plot range
        mask = np.logical_and(alllinelocs > xlim[0], alllinelocs < xlim[1])
        if sum(mask) != 0:
            for loc in linelocs:
                #   Align parameter
                align = loc[1]
                #   Wave length
                wave = loc[0]

                if not xlim[0] < wave < xlim[1]:
                    continue

                #   Plot identifier
                trans = axs.get_xaxis_transform()
                ann = axs.annotate(
                    linestr,
                    xy=(wave, 1.05),
                    xycoords=trans,
                    ha=align,
                )
                axs.axvline(
                    wave,
                    color="lightgrey",
                    linestyle="--",
                    zorder=1,
                )


def plot_merged(wave_merged, flux_merged, obj_name, normalize, mtype=''):
    """
        Plot merged spectrum

        Parameters
        ----------
        wave_merged             : `numpy.ndarray`
            Wavelength data

        flux_merged             : `numpy.ndarray`
            Flux data

        obj_name                : `string`
            Name of the object

        normalize               : `boolean`
            If `True`, it is assumed that the flux is normalized.

        mtype                   : `string`, optional
            String that characterize the order merging procedure. It that will
            be part of the file names.
            Default is ``''``.

    """
    print(
        f'      Plot total spectrum:{Bcolors.OKBLUE} {obj_name}{Bcolors.ENDC}'
    )

    #   Define figure
    fig1 = plt.figure(figsize=(10, 5))

    #   Set plot range
    offset = (max(flux_merged) - min(flux_merged)) * 0.02
    plt.ylim([min(flux_merged) - offset, max(flux_merged) + offset])
    offset = (max(wave_merged) - min(wave_merged)) * 0.02
    plt.xlim([min(wave_merged) - offset, max(wave_merged) + offset])

    #   Set title and label etc.
    plt.suptitle(f"{obj_name.replace('_', ' ')} (orders merged by {mtype})")
    plt.xlabel(r'Wavelength $\lambda\,[\AA]$')
    if normalize:
        plt.ylabel('Normalized flux')

        #   Set plot range
        minflux = np.max([np.min(flux_merged) * 0.94, 0])
        maxflux = np.min([np.max(flux_merged) * 1.06, 2])
        plt.ylim([minflux, maxflux])
    else:
        plt.ylabel('Relative flux')
        minflux = np.max([np.min(flux_merged) * 0.94, 0])
        maxflux = np.max(flux_merged) * 1.06
        plt.ylim([minflux, maxflux])

    plt.tick_params(top=True, right=True, which='both', direction='in')
    plt.minorticks_on()

    #   Plot spectrum
    plt.plot(wave_merged, flux_merged, 'b-', linewidth=0.5)

    #   Save spectrum as PDF file
    print(
        f'      Create spectrum plot{Bcolors.OKBLUE} output/'
        f'spectrum_total_{mtype}-merged_{obj_name}.pdf {Bcolors.ENDC}'
    )
    os.system("mkdir -p output")
    pp = PdfPages(
        f'output/spectrum_total_{mtype}-merged_{obj_name}.pdf'
    )

    pp.savefig(fig1, dpi=300, transparent=True)
    pp.close()
    plt.clf()
    plt.close()

    #   Write data to CSV file
    print(
        f'      Write data to{Bcolors.OKBLUE} output/{obj_name}_{mtype}'
        f'-merged_spectrum_total.csv {Bcolors.ENDC}'
    )

    np.savetxt(
        f"output/{obj_name}_{mtype}-merged_spectrum_total.csv",
        np.transpose((wave_merged, flux_merged)),
        delimiter=",",
    )
    np.savetxt(
        f"output/{obj_name}_{mtype}-merged_spectrum_total.dat",
        np.transpose((wave_merged, flux_merged)),
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


def rename_elem(string):
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
        return (string_split[0] + to_roman(string_split[1])).strip()
    else:
        return string


def generate_lines_from_file(wl, flx, ions, lineFile="atomic_lines.tsv",
                             line_lower_continuum=3.):
    """
        Read line file and prepare lines of the ions that are requested.
        Gets the lines with the larges gaunt factors. The lines need to be
        within the plot range.

        Parameters
        ----------
        wl                      : `numpy.ndarray`
            Wavelength data

        flx                     : `numpy.ndarray`
            Flux data

        ions                    : `list` of `string`
            Ions whose lines are to be extracted from the file.

        lineFile                : `string`, optional
            Name/Path to the file with spectral line data
            Default is ``atomic_lines.tsv``.

        line_lower_continuum    : `float`, optional
            Percent the line flux must be lower than the continuum
            Default is ``3``.

        Returns
        -------
        lines               : `dict`
            Line information: key=ion identifier, value=wavelength
    """
    #   Interpolate on wavelength and flux to allow determination of flux
    #   level at the individual line positions
    f = interp1d(wl, flx, fill_value=1, bounds_error=False)

    #   Read line file and restrict data to a reasonable range
    try:
        linefile = pd.read_csv(lineFile, delimiter="\t")
    except:
        print(
            f"{Bcolors.WARNING}   Line file not found. Check variable "
            f"'lineFile'. Specified was {lineFile}. {Bcolors.ENDC}"
        )
        return {}
    linefile = linefile.loc[
        (linefile["wave_A"] > 4200) & (linefile["wave_A"] < 7600)
        ]

    #   Replace Arabic numbers with Roman numbers
    linefile['element'] = linefile['element'].apply(rename_elem)

    #   Setup lists and a dict for the line data
    lines = {}
    alllines = []
    allstrengths = []
    allions = []

    #   Restrict line data to strong lines and lines with the largest
    #   gaunt factors
    for ion in ions:
        #   Restrict line data to current ion and the 25 strongest lines
        subset = linefile.loc[linefile["element"] == ion]
        subset = subset.sort_values('loggf').tail(25)

        for ind, row in subset.iterrows():
            #   Get line wavelength and index of closest wavelength in wl array
            linewl = float(row["wave_A"])
            index = np.argmin(np.abs(wl - linewl))

            #   Mean flux around line wavelength
            try:
                mean_flx = np.mean(flx[index - 40:index + 40])
            except IndexError:
                continue

            #   Skip weak (not deep) lines
            if f(linewl) > (100 - line_lower_continuum) / 100 * mean_flx:
                if row["element"] != 'HI':
                    continue

            #   Remove/skip lines with the same wavelength if there is one
            #   with a larger gaunt factor
            if linewl in alllines:
                ind = alllines.index(linewl)
                if row["loggf"] > allstrengths[ind]:
                    lines[allions[ind]].remove([linewl, "center"])

                    alllines.pop(ind)
                    allstrengths.pop(ind)
                    allions.pop(ind)
                else:
                    continue

            #   Check if ion already exists in the `lines` dictionary
            try:
                lines[row["element"]]
            except KeyError:
                lines[row["element"]] = []

            #   Fill list with line data
            lines[row["element"]].append([linewl, "center"])
            alllines.append(linewl)
            allstrengths.append(float(row["loggf"]))
            allions.append(row["element"])

    #   Convert lists to arrays
    alllines = np.array(alllines)
    allstrengths = np.array(allstrengths)
    allions = np.array(allions)

    #   Remove lines that are too close to each other - keep strongest
    for k, line in enumerate(alllines):
        #   Find close lines
        mask = np.logical_and(line - 4 < alllines, alllines < line + 4)
        if sum(mask) > 1:
            #   Hydrogen is always the strongest line
            if "HI" in allions[np.argwhere(alllines == line)]:
                strongest = line
            else:
                #   Find strong line based on gaunt factor
                index_strongest = np.argmax(allstrengths[mask])
                strongest = alllines[mask][index_strongest]

            #   Get other lines (weaker lines)
            otherlines = alllines[mask]
            otherlines = otherlines[otherlines != strongest]

            #   Remove other lines
            for otherline in otherlines:
                i = np.argwhere(alllines == otherline)[0][0]
                lines[allions[i]].remove([alllines[i], "center"])

                alllines = np.delete(alllines, i)
                allstrengths = np.delete(allstrengths, i)
                allions = np.delete(allions, i)

    #   Set horizontal alignment for identifier string, if lines are to close
    for ion, llist in lines.items():
        for l in llist:
            #   Line wave length
            wave = l[0]

            #   Find close lines
            # if sum(np.logical_and(alllines > l - 0.5, alllines < l + 0.5)) > 1:
            # np.delete(alllines, alllines == l)
            # lines[element].remove([l, "center"])
            # continue

            #   Identify close lines
            mask = np.logical_and(alllines > wave - 10, alllines < wave + 10)
            if sum(mask) > 1:
                otherlines = alllines[mask]
                otherlines = otherlines[otherlines != wave]

                #   Set alignment according to line wavelength
                for otherline in otherlines:
                    index = np.argwhere(alllines == otherline)[0][0]
                    other_ion = allions[index]
                    other_ion_waves = np.array(lines[other_ion])[:, 0]
                    index_wave = list(other_ion_waves).index(str(otherline))
                    if otherline > wave:
                        lines[other_ion][index_wave][1] = "left"
                    else:
                        lines[other_ion][index_wave][1] = "right"

    return lines


def plot_panels(wave_merged, flux_merged, obj_name, normalize, lineFile,
                panel_wave_range, mtype=''):
    """
        Plot merged data in individual panels and create PDFs

        Parameters
        ----------
        wave_merged             : `numpy.ndarray`
            Wavelength data

        flux_merged             : `numpy.ndarray`
            Flux data

        obj_name              : `string`
            Name of the object

        normalize               : `boolean`
            If `True`, it is assumed that the flux is normalized.

        lineFile                : `string`
            Path to the file with the line identifications

        panel_wave_range        : `float` or `integer`
            Wavelength range for the individual panels

        mtype                   : `string`, optional
            String that characterize the order merging procedure. It that will
            be part of the file names.
            Default is ``''``.
    """
    print(
        f"      Plot individual panels:{Bcolors.OKBLUE} {obj_name}{Bcolors.ENDC}"
    )

    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    N = len(wave_merged)
    i = 1
    j_p = 0
    for j in range(0, N):
        if wave_merged[j] >= wave_merged[0] + panel_wave_range * i or j + 1 == N:

            #   Set wave and flux range
            wave_range = wave_merged[j_p:j]
            flux_range = flux_merged[j_p:j]

            #   Define figure and set labels etc.
            fig = plt.figure(figsize=(12, 6))
            plt.suptitle(
                '{} ({} - {}$\,\AA$, orders merged by {})'.format(
                    obj_name.replace('_', ' '),
                    int(wave_range[0]),
                    int(wave_range[-1]),
                    mtype,
                )
            )
            plt.xlabel(r'Wavelength $\lambda\,[\AA]$')
            if normalize:
                plt.ylabel('Normalized flux')

                #   Set plot range
                minflux = np.max([np.min(flux_range) * 0.94, 0])
                maxflux = np.min([np.max(flux_range) * 1.06, 2])
                plt.ylim([minflux, maxflux])
            else:
                plt.ylabel('Relative flux')
                minflux = np.max([np.min(flux_range) * 0.94, 0])
                maxflux = np.max(flux_range) * 1.06
                plt.ylim([minflux, maxflux])

            plt.tick_params(top=True, right=True, which='both', direction='in')
            plt.minorticks_on()

            #   Plot idents
            add_idents(wave_range, flux_range, lineFile, minflux, maxflux)

            #   Plot spectrum
            plt.step(wave_range, flux_range, 'b-', linewidth=0.5)

            #   Save panel as PDF file
            plt.savefig(
                f'{temp_dir.name}/{obj_name}_{mtype}-merged_{i}.pdf',
                bbox_inches='tight',
            )
            plt.close()

            #   Write data to CSV file
            print(
                f'      Write data to{Bcolors.OKBLUE} output/{obj_name}_'
                f'{mtype}-merged_spectrum_{i}.csv {Bcolors.ENDC}'
            )
            np.savetxt(
                f"output/{obj_name}_{mtype}-merged_spectrum_{i}.csv",
                np.transpose((wave_range, flux_range)),
                delimiter=",",
            )

            i += 1
            j_p = copy.deepcopy(j)

    #   Merge individual plotted bins with pdfunite
    print(
        f"      Create{Bcolors.OKBLUE} output/spectrum_panels_{mtype}-merged_"
        f"{obj_name}.pdf {Bcolors.ENDC}"
    )

    merger = PdfMerger()
    file_path = f"output/spectrum_panels_{mtype}-merged_{obj_name}.pdf"
    for f in os.listdir(temp_dir.name):
        merger.append(f"{temp_dir.name}/{f}")
    with open(file_path, "wb") as new_file:
        merger.write(new_file)

    merger.close()


def plot_panels_fabian(wave, flux, obj_name, normalize, lines, divinto=21,
                       mtype='', rv=0., cont=None, n=5):
    """
        Plot merged data in individual panels and create PDFs

        Parameters
        ----------
        wave                     : `numpy.ndarray`
            Wavelength data

        flux                    : `numpy.ndarray`
            Flux data

        obj_name                : `string`
            Name of the object

        normalize               : `boolean`
            If `True`, it is assumed that the flux is normalized.

        lines                   : `dictionary`
            Line information: key=line identifier, value[0]=wavelength,
                                                   value[1]=alignment parameter

        divinto                 : `integer`, optional
            Number of panels that should be used to display the spectrum
            Default is ``21``.

        mtype                   : `string`, optional
            String that characterize the order merging procedure. It that will
            be part of the file names.
            Default is ``''``.

        rv                      : `float`, optional
            Radial velocity

        cont                    : `numpy.ndarray`
            Flux data of the continuum

        n                       : `integer`, optional
            Number of panels on each page/plot.
            Default is ``5``.
    """
    print(
        f"      Plot individual panels:{Bcolors.OKBLUE} {obj_name}{Bcolors.ENDC}"
    )

    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Limit to a range which can be divide by ``divinto``
    while len(flux) % divinto != 0:
        wave = wave[:-1]
        flux = flux[:-1]
        if cont is not None:
            cont = cont[:-1]

    #   Split data into ``divinto`` ranges
    wave_s = np.split(wave, divinto)
    flux_s = np.split(flux, divinto)
    if cont is not None:
        cont_s = np.split(cont, divinto)
    else:
        cont_s = [0 for l in wave_s]

    #   Initialize plot
    fig, axs = plt.subplots(n, 1, figsize=(8.2, 11.6))

    for i, [w, f, c] in enumerate(zip(wave_s, flux_s, cont_s)):
        #   Axis index
        axind = i
        while axind >= 0:
            axind -= n
        i += 1

        # plt.tick_params(top=True, right=True, which='both', direction='in')
        # plt.minorticks_on()

        #   Plot data
        axs[axind].step(w, f, color='#0066ff', zorder=5)

        #   Plot continuum
        if cont is not None:
            axs[axind].plot(w, c, color='red', zorder=5)

        #   Set ranges
        xlim = (np.amin(w), np.amax(w))
        axs[axind].set_xlim(xlim)
        if normalize:
            minflux = np.max([np.amin(f) * 0.94, 0])
            minflux = 0.
            maxflux = np.min([np.amax(f) * 1.06, 2])
            axs[axind].set_ylim([minflux, maxflux])
            axs[axind].set_ylabel('Normalized flux')
            axs[axind].legend(["Measured"], loc="lower right")
            axs[axind].axline(
                (xlim[0], 1.0),
                xy2=(xlim[1], 1.0),
                color="darkgrey",
                linestyle="-",
                linewidth=0.5,
                zorder=1,
            )
        else:
            minflux = np.max([np.amin(f) * 0.94, 0])
            maxflux = np.amax(f) * 1.06
            axs[axind].set_ylim([minflux, maxflux])
            axs[axind].set_ylabel('Relative flux')
            if cont is not None:
                axs[axind].legend(["Measured", "Continuum"], loc="lower right")
            else:
                axs[axind].legend(["Measured"], loc="lower right")

        #   Set plot title, labels, etc.
        fig.suptitle(
            f"Spectrum of {obj_name} ({mtype}) - Page "
            f"[{math.ceil(i / n)}/{math.ceil(divinto / n)}]",
            fontsize=15,
        )

        axs[n - 1].set_xlabel("$\lambda$ [Å]")

        if ((i - 1) % n) == 0:
            axs[axind].annotate(
                f"Corrected for radial velocity "
                f"$v_{{\mathrm{{rad}}}}={round(rv, 2)}$ km/s",
                xy=(axs[axind].get_xlim()[0], 1.125),
                xycoords=axs[axind].get_xaxis_transform(),
                ha="left",
                fontsize=10,
            )

        #   Plot line identifier
        add_idents_fabian(lines, xlim, axs[axind])

        #   Save plot
        if i % n == 0:
            plt.tight_layout()
            plt.savefig(f"{temp_dir.name}/{int(i / n)}.pdf")
            plt.close()
            fig, axs = plt.subplots(n, 1, figsize=(8.2, 11.6))

            ##   Write data to CSV file
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
        f"      Create{Bcolors.OKBLUE} output/spectrum_panels_{mtype}-merged_"
        f"{obj_name}.pdf {Bcolors.ENDC}"
    )
    merger = PdfMerger()
    file_path = f"output/spectrum_panels_{mtype}-merged_{obj_name}.pdf"
    for f in os.listdir(temp_dir.name):
        merger.append(f"{temp_dir.name}/{f}")
    with open(file_path, "wb") as new_file:
        merger.write(new_file)

    merger.close()


############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    ###
    #   Prepare stuff
    #
    #  Sanitize object name
    objectname = objectname.replace(' ', '_')

    #  Create output directory
    os.system("mkdir -p output")

    ########################################################################
    ####                      MIDAS merged spectrum                      ###
    ########################################################################
    print(
        f"{Bcolors.BOLD}Processing spectrum merged by MIDAS{Bcolors.ENDC}"
    )

    ###
    #   Load spectra
    #
    #   Select FIT or fit
    File_merged = check_file_name(File_merged)

    #   Open fits file and read the data section and header section into a
    #   2 dimensional array
    print(f"{Bcolors.BOLD}   Read files ...{Bcolors.ENDC}")
    wave_merged, flux_merged = read_baches_merged(File_merged)

    ###
    #   Normalize flux
    #
    if normalize:
        print(
            f"{Bcolors.BOLD}   Normalize merged MIDAS spectrum{Bcolors.ENDC}"
        )
        if norm_version == 'specutils_continuum':
            wave_merged, flux_merged = norm_spectrum_interval(
                wave_merged,
                flux_merged,
                median_window,
                porder,
            )
            wave_merged = wave_merged.value
            flux_merged = flux_merged.value
            cont = None
        elif norm_version == 'median_max_window':
            wave_merged = wave_merged.value
            flux_merged = flux_merged.value
            flux_merged, cont = normalize_spectrum_fabian(
                wave_merged,
                flux_merged,
                name=objectname,
            )
        else:
            print(
                f"{Bcolors.FAIL}   Normalize method not known. Check variable:"
                f" 'norm_version'. Skipping normalization!{Bcolors.ENDC}"
            )
            cont = None
    else:
        wave_merged = wave_merged.value
        flux_merged = flux_merged.value
        cont = None

    ###
    #   Correct for radial velocity -> 1) calculate barycentric velocity
    #                                     correction
    #                                  2) add radial velocity
    #
    velocity_correction = radial_velocity
    if correct_bary:
        try:
            bvc = bary_correction(File_merged)
        except RuntimeError:
            print(
                f"{Bcolors.WARNING}   Barycentric velocity correction could not "
                f"be determined. Assume 0 km/s.{Bcolors.ENDC}"
            )
            bvc = 0.
        velocity_correction = radial_velocity - bvc

    wave_merged = correct_for_rv(
        wave_merged,
        flux_merged,
        velocity_correction,
        objectname,
        # plot=True,
        mtype='MIDAS',
    )

    ###
    #   Plot merged data in one plot
    #
    print(f"{Bcolors.BOLD}   Plot spectra{Bcolors.ENDC}")
    plot_merged(
        wave_merged,
        flux_merged,
        objectname,
        normalize,
        mtype='MIDAS',
    )

    ###
    #   Plot merged data in individual panels
    #
    if panel_version == 'old':
        plot_panels(
            wave_merged,
            flux_merged,
            objectname,
            normalize,
            lineFile,
            panel_wave_range,
            mtype='MIDAS',
        )
    elif panel_version == 'default':
        #   Generate ident lines
        lines = generate_lines_from_file(
            wave_merged,
            flux_merged,
            ions,
            lineFile=lineFile,
            line_lower_continuum=line_lower_continuum,
        )

        plot_panels_fabian(
            wave_merged,
            flux_merged,
            objectname,
            normalize,
            {**lines, **manual_lines},
            divinto=divinto,
            mtype='MIDAS',
            rv=radial_velocity,
            cont=cont,
            n=n_panel,
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
        File_orders = check_file_name(File_orders)

        #   Open fits file and read the data and header section
        print(f"{Bcolors.BOLD}   Read files ...{Bcolors.ENDC}.")
        wave_list, flux_list = read_baches(File_orders)

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
                orders.append(Spectrum1D(spectral_axis=wave, flux=flux_list[i]))

            #   Normalize & merge spectra
            merged_wave, merged_flux = norm_merge_spectra(
                orders,
                median_window=median_window,
                order=porder,
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
                trim_value=trim_value,
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
                merged_wave, merged_flux = norm_spectrum_interval(
                    merged_wave,
                    merged_flux,
                    median_window,
                    porder,
                )
                merged_wave = merged_wave.value
                merged_flux = merged_flux.value
                cont = None
            elif norm_version == 'median_max_window':
                merged_wave = merged_wave.value
                merged_flux = merged_flux.value
                merged_flux, cont = normalize_spectrum_fabian(
                    merged_wave,
                    merged_flux,
                    name=objectname,
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
            cont = None

        ###
        #   Correct for radial velocity
        #
        merged_wave = correct_for_rv(
            merged_wave,
            merged_flux,
            velocity_correction,
            objectname,
            # plot=True,
            mtype='PYTHON',
        )

        ###
        #   Plot merged data in one plot
        #
        print(f"{Bcolors.BOLD}   Plot spectra{Bcolors.ENDC}")
        plot_merged(
            merged_wave,
            merged_flux,
            objectname,
            normalize,
            mtype='PYTHON',
        )

        ###
        #   Plot merged data in individual panels
        #
        if panel_version == 'old':
            plot_panels(
                merged_wave,
                merged_flux,
                objectname,
                normalize,
                lineFile,
                panel_wave_range,
                mtype='PYTHON',
            )
        elif panel_version == 'default':
            #   Generate ident lines
            lines = generate_lines_from_file(
                merged_wave,
                merged_flux,
                ions,
                lineFile=lineFile,
                line_lower_continuum=line_lower_continuum,
            )

            plot_panels_fabian(
                merged_wave,
                merged_flux,
                objectname,
                normalize,
                {**lines, **manual_lines},
                divinto=divinto,
                mtype='PYTHON',
                rv=radial_velocity,
                cont=cont,
                n=n_panel,
            )
        else:
            print(
                f"{Bcolors.FAIL}   Version for panel version not known. "
                f"Check variable: 'panel_version'. Skipping panel "
                f"plot!{Bcolors.ENDC}"
            )

    print(f'{Bcolors.OKGREEN}DONE{Bcolors.ENDC}')
