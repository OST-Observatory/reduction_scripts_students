# ! /usr/bin/python3
# -*- coding: utf-8 -*-

##################################################################################
#                               Script Parameters                                #
##################################################################################

#   Name of the object
object_name: str = "T_CrB"

#   Radial velocity [km/s]
#       The specification of the radial velocity is necessary for the
#       line identification to work correctly (see below).
radial_velocity: float = 0.


###
#   Extraction regions
#
#   Region containing the science spectrum
spec_region_start: int = 759
spec_region_end: int = 772

#   Sky background region (inside the slit)
background_sky_start: int = 710
background_sky_end: int = 730


###
#   Plot range
#
#   Set the variables to '?' for an automatic resizing
lambda_min: str | float = '?'
lambda_max: str | float = '?'
lambda_min: str | float = 4000.
lambda_max: str | float = 8000.
# lambda_max: str | float = 9500.


###
#   Normalization ?
#   Possibilities: True or False
#
normalize: bool = False
# normalize: bool = True


###
#   Line identifications
#
#   Ions for which line markers are to be drawn.
#   Example: ["HI", "FeI", ...]
ions: list[str] = ["HI",]

#   Add lines that ar not in the line file
#   Format: {"Element descriptor": [[wavelength, alignment parameter]]}
#           alignment parameter possibilities: "center", "left", "right"
manual_lines: dict[str, list[list[float | str]]] = {"Example Element": [[0., "center"]]}

#   Percent the line flux must be lower than the continuum
percentage_line_flux_must_be_below_continuum: float = 3.


##################################################################################
#          The following parameters usually do not need to be adjusted           #
##################################################################################

###
#   Files
#
#   Science spectrum file
spectrum_file: str = 'master_spectrum.fit'

#   Flat spectrum file
flat_file: str = 'master_flat.fit'

#   Calibration file
calibration_file: str = 'calibration_spectrum.dat'


###
#   Apply Flat calibration
#
apply_flat_calibration: bool = True
# apply_flat_calibration: bool = False


###
#   Normalization options
#
#   Normalization version
#   Possibilities: `median_max_window` or `specutils_continuum`
#   Default is `median_max_window`
# norm_version: str = 'specutils_continuum'
norm_version: str = 'median_max_window'

#   Normalize before merging of the orders
#   Possibilities: True or False
#   If `False` the orders will be merged first and the merged spectrum will
#   be normalized afterward.
#   (norm_version = specutils_continuum)
norm_individual: bool = False
# norm_individual: bool = True

#   Order of the polynomial used to normalize the spectra
#   (norm_version = specutils_continuum)
porder: int = 5

#   Wavelength window used to estimate continuum points
#   (norm_version = specutils_continuum)
median_window: int = 61

#   Number of wavelength bins that should be removed from the beginning
#   of the orders to improve flux normalization between the orders
trim_value_start_order: int = 0
# trim_value_start_order: int = 400


###
#   Panel plot options
#
#   Number of panels on each page/plot
n_panels_per_page: int = 4

#   Number of panels into which the spectrum should be split
n_panels: int = 4


###
#   Line identifications
#
#   File containing line identifications
# line_file: str = ""
# line_file: str = "absorption_lines.dat"
line_file: str = "~/projects/reduction_scripts_students/n1_baches/atomic_lines.tsv"


###
#   Apply barycentric correction?
#
barycenter_correction: bool = False
barycenter_correction: bool = True


###
# Data output
#
spectrum_output_file: str = 'stern_spectrum.pdf'
spectrum_output_data_file: str = 'stern_spectrum.dat'


##################################################################################
#                            Script Routines                                     #
##################################################################################

import os

import tempfile

import math

import numpy as np

from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
from scipy.ndimage import median_filter, maximum_filter
from scipy import interpolate
from scipy.constants import c

import pandas as pd

from astropy.io import fits
from astropy.nddata import CCDData
from astropy.stats import sigma_clip, mad_std
from astropy.modeling import models, fitting
from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time
from astropy import uncertainty as unc
import astropy.units as u

from specutils import Spectrum1D
from specutils.fitting import fit_generic_continuum
from specutils.spectra import SpectralRegion

from PyPDF2 import PdfMerger

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42


class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


##################################################################################
#                               Functions                                        #
##################################################################################

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


def extract_spectrum_data(
        image_data: CCDData, row_start: int | str, row_end: int | str) -> None:
    """
    Extracts the spectrum data from a number of rows specified by 'row_start'
    and 'row_end'.

    Parameter
    ---------
    image_data
        2D image data containing the spectral information

    row_start
        Row from which the extraction will start

    row_end
        Row where the extraction will end

    Returns:
    --------
    spectrum_1d
        Distribution with 1D flux data of the spectrum
    """
    #   Get data from image in requested row range
    spectrum_ccd = image_data[row_start:row_end, :]

    #   Set up distribution with the data
    distribution_samples: int = 1000
    spectrum_2d = unc.normal(
        spectrum_ccd.data,
        std=spectrum_ccd.uncertainty.array,
        n_samples=distribution_samples,
    )

    #   Calculate the median over the range of rows
    spectrum_1d = np.median(spectrum_2d, axis=0)

    # spectrum_mask_1d = np.sum(spectrum_ccd.mask, axis=0)
    # print(type(spectrum_mask_1d))
    # print(np.any(spectrum_mask_1d))
    # print(np.any(spectrum_ccd.mask))

    return spectrum_1d


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
        if merge_type:
            tile_first_part = f"Spectrum of {object_name} ({merge_type}) - Page "
        else:
            tile_first_part = f"Spectrum of {object_name} - Page "
        figure.suptitle(
            f"{tile_first_part}"
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
    if merge_type:
        print(
            f"      Create{Bcolors.OKBLUE} output/spectrum_panels_{merge_type}-"
            f"merged_{object_name}.pdf {Bcolors.ENDC}"
        )

        file_path = f"output/spectrum_panels_{merge_type}-merged_{object_name}.pdf"
    else:
        print(
            f"      Create{Bcolors.OKBLUE} output/spectrum_panels_"
            f"{object_name}.pdf {Bcolors.ENDC}"
        )

        file_path = f"output/spectrum_panels_{object_name}.pdf"

    merger = PdfMerger()
    for flux_element in os.listdir(temp_dir.name):
        merger.append(f"{temp_dir.name}/{flux_element}")
    with open(file_path, "wb") as new_file:
        merger.write(new_file)

    merger.close()


def plot_merged(wave_length_merged, flux_merged, object_name, normalize,
                lines_to_mark, merge_type=''):
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

        lines_to_mark           : `dictionary`
            Line information: key=line identifier, value[0]=wavelength,
                                                   value[1]=alignment parameter

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
    # figure = plt.figure(figsize=(10, 5))
    # figure = plt.figure(figsize=(11.69, 8.27))
    figure = plt.figure(figsize=(11.69, 6.))

    #   Set plot range
    flux_offset = (max(flux_merged) - min(flux_merged)) * 0.02
    plt.ylim([min(flux_merged) - flux_offset, max(flux_merged) + flux_offset])
    wave_length_offset = (max(wave_length_merged) - min(wave_length_merged)) * 0.02
    plt.xlim([
        min(wave_length_merged) - wave_length_offset,
        max(wave_length_merged) + wave_length_offset
    ])

    #   Set title and label etc.
    if merge_type:
        tile_part_two = f"(orders merged by {merge_type})"
    else:
        tile_part_two = ""
    plt.suptitle(
        f"{object_name.replace('_', ' ')} {tile_part_two}"
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

    #   Plot line identifier
    add_idents_fabian(
        lines_to_mark,
        (np.amin(wave_length_merged), np.amax(wave_length_merged)),
        plt.gca(),
    )

    #   Save spectrum as PDF file
    os.system("mkdir -p output")
    if merge_type:
        print(
            f'      Create spectrum plot{Bcolors.OKBLUE} output/'
            f'spectrum_total_{merge_type}-merged_{object_name}.pdf {Bcolors.ENDC}'
        )
        pp = PdfPages(
            f'output/spectrum_total_{merge_type}-merged_{object_name}.pdf'
        )
    else:
        print(
            f'      Create spectrum plot{Bcolors.OKBLUE} output/'
            f'spectrum_total_{object_name}.pdf {Bcolors.ENDC}'
        )
        pp = PdfPages(
            f'output/spectrum_total_{object_name}.pdf'
        )

    pp.savefig(figure, dpi=300, transparent=True)
    pp.close()
    plt.clf()
    plt.close()

    #   Write data to CSV file
    if merge_type:
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
    else:
        print(
            f'      Write data to{Bcolors.OKBLUE} output/{object_name}_'
            f'spectrum_total.csv {Bcolors.ENDC}'
        )

        np.savetxt(
            f"output/{object_name}_spectrum_total.csv",
            np.transpose((wave_length_merged, flux_merged)),
            delimiter=",",
        )
        np.savetxt(
            f"output/{object_name}_spectrum_total.dat",
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

            #   Get other lines next to the current line (weaker lines)
            other_lines_wave_length = line_wave_length_array[mask]
            other_lines_wave_length = other_lines_wave_length[other_lines_wave_length != strongest_line]

            #   Remove the other lines
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
                line_strength_array = np.delete(
                    line_strength_array,
                    index_other_line
                )
                # line_strength_list = np.delete(line_strength_array, index_other_line)
                ion_array = np.delete(ion_array, index_other_line)

    #   Set horizontal alignment for identifier string, if lines are too close
    for ion, line_list in line_dict.items():
        for line in line_list:
            #   Line wave length
            line_wave_length = line[0]

            #   Find close lines
            # if sum(np.logical_and(line_wave_length_array > l - 0.5, line_wave_length_array < l + 0.5)) > 1:
            # np.delete(line_wave_length_array, line_wave_length_array == l)
            # line_dict[element].remove([l, "center"])
            # continue

            #   Identify close lines
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

    medium_window_size = 50
    max_window_size = 400
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


def main() -> None:
    """
    Main function performing data extraction and calibration

    """
    print(Bcolors.BOLD + "   Extract spectral data" + Bcolors.ENDC)
    #   Load 2D spectras
    ccd_spectrum = CCDData.read(spectrum_file)
    ccd_flat = CCDData.read(flat_file)

    #   Extract the spectrum data
    object_flux = extract_spectrum_data(
        ccd_spectrum,
        spec_region_start,
        spec_region_end,
    )

    #   Extract sky spectrum
    sky_flux = extract_spectrum_data(
        ccd_spectrum,
        background_sky_start,
        background_sky_end,
    )

    #   Extract flat spectrum
    flat_flux = extract_spectrum_data(
        ccd_flat,
        spec_region_start,
        spec_region_end,
    )

    #   Remove sky spectrum from object spectrum
    flux = object_flux - sky_flux

    #   Normalize and apply flat
    if apply_flat_calibration:
        normalized_flat = flat_flux / np.max(flat_flux)
        flux /= normalized_flat

    #   Get data from distribution
    flux = flux.pdf_median()

    print(Bcolors.BOLD + "   Perform wavelength calibration" + Bcolors.ENDC)
    calibration_data_file = open(calibration_file, 'r')
    wavelength = []
    for line in calibration_data_file:
        liste = line.split()
        if len(liste) == 0:
            continue
        if lambda_min != '?' and float(liste[0]) < lambda_min:
            flux = flux[1:]
            continue
        if lambda_max != '?' and float(liste[0]) > lambda_max:
            flux = flux[:-1]
            continue
        wavelength.append(float(liste[0]))
    wavelength = np.asarray(wavelength)

    print(Bcolors.BOLD + "   Create spectral plot " + Bcolors.OKBLUE + spectrum_output_file + Bcolors.ENDC)

    #   Normalize flux
    if normalize:
        print(
            f"{Bcolors.BOLD}   Normalize merged spectrum{Bcolors.ENDC}"
        )
        if norm_version == 'specutils_continuum':
            wavelength, flux = normalize_spectrum_interval(
                wavelength,
                flux,
                median_window,
                porder,
            )
            wavelength = wavelength.value
            flux = flux.value
            continuum = None
        elif norm_version == 'median_max_window':
            flux, continuum = normalize_spectrum_fabian(
                wavelength,
                flux,
                object_name=object_name,
            )
        else:
            print(
                f"{Bcolors.FAIL}   Normalize method not known. Check variable:"
                f" 'norm_version'. Skipping normalization!{Bcolors.ENDC}"
            )
            continuum = None
    else:
        continuum = None

    ###
    #   Correct for radial velocity -> 1) calculate barycentric velocity
    #                                     correction
    #                                  2) add radial velocity
    #
    velocity_correction = radial_velocity
    if barycenter_correction:
        try:
            bvc = bary_correction(spectrum_file)
        except RuntimeError:
            print(
                f"{Bcolors.WARNING}   Barycentric velocity correction could "
                f"not be determined. Assume 0 km/s.{Bcolors.ENDC}"
            )
            bvc = 0.
        velocity_correction = radial_velocity - bvc

    wavelength = correct_for_radial_velocity(
        wavelength,
        flux,
        velocity_correction,
        object_name,
        # plot=True,
    )

    ###
    #   Generate ident lines
    #
    lines = generate_lines_from_file(
        wavelength,
        flux,
        ions,
        line_file=line_file,
        percentage_line_flux_must_be_below_continuum=percentage_line_flux_must_be_below_continuum,
    )

    ###
    #   Plot merged data in one plot
    #
    print(f"{Bcolors.BOLD}   Plot spectra{Bcolors.ENDC}")
    plot_merged(
        wavelength,
        flux,
        object_name,
        normalize,
        {**lines, **manual_lines},
    )

    ###
    #   Plot merged data in individual panels
    #
    plot_panels_fabian(
        wavelength,
        flux,
        object_name,
        normalize,
        {**lines, **manual_lines},
        n_panels=n_panels,
        radial_velocity=radial_velocity,
        continuum=continuum,
        n_panels_per_page=n_panels_per_page,
    )


############################################################################
#                                  Main                                    #
############################################################################

if __name__ == '__main__':
    main()

#
# # spectrum star
# spec = []
# # spectrum flatfield
# spec2 = []
# # Sky background spectrum
# skyspec = []
#
# for i in range(0, len(sience_reduced_array[0])):
#     spec_tmp, bg_tmp = 0, 0
#     for j in range(spec_region_start, spec_region_end):
#         spec_tmp += sience_reduced_array[j][i]
#     # for k in range(bgRegionStart,bgRegionEnd):
#     # bg_tmp += sience_reduced_array[k][i]
#     spec_tmp /= abs((spec_region_end - spec_region_start))
#     # bg_tmp /=  abs((bgRegionEnd-bgRegionStart))
#     spec.append(spec_tmp - bg_tmp)

# for i in range(0, len(flatfield_out_array[0])):
#     spec_tmp, bg_tmp = 0, 0
#     for j in range(spec_region_start, spec_region_end):
#         spec_tmp += flatfield_out_array[j][i]
#     # for k in range(bgRegionStart,bgRegionEnd):
#     # bg_tmp += flatfield_out_array[k][i]
#     spec_tmp /= abs((spec_region_end - spec_region_start))
#     # bg_tmp /=  abs((bgRegionEnd-bgRegionStart))
#     spec2.append(spec_tmp - bg_tmp)
#
# for i in range(0, len(sience_reduced_array[0])):
#     skybg, bg_tmp = 0, 0
#     for j in range(background_sky_start, background_sky_end):
#         skybg += sience_reduced_array[j][i]
#     # for k in range(bgRegionStart,bgRegionEnd):
#     # bg_tmp += sience_reduced_array[k][i]
#     skybg /= abs((background_sky_start - background_sky_end))
#     # bg_tmp /=  abs((bgRegionEnd-bgRegionStart))
#     skyspec.append(skybg - bg_tmp)
#
# print(bcolors.BOLD + "   Apply Flatfield correction to science spectrum" + bcolors.ENDC)
#
# # normalize the flatfield
# spec2 = np.asarray(spec2) / max(spec2)
# # apply flatfield to science spectrum and sky spectrum
# spec = spec / spec2
# skyspec = skyspec / spec2
#
# print(bcolors.BOLD + "   Perform wavelength calibration" + bcolors.ENDC)
#
# calibfile = open(calibration_file, 'r')
# wavelengthrange = []
# for line in calibfile:
#     liste = line.split()
#     if len(liste) == 0:
#         continue
#     if lambda_min != '?' and float(liste[0]) < lambda_min:
#         spec = spec[1:]
#         spec2 = spec2[1:]
#         skyspec = skyspec[1:]
#         continue
#     if lambda_max != '?' and float(liste[0]) > lambda_max:
#         spec = spec[:-1]
#         spec2 = spec2[:-1]
#         skyspec = skyspec[:-1]
#         continue
#     wavelengthrange.append(float(liste[0]))
#
# print(bcolors.BOLD + "   Create spectral plot " + bcolors.OKBLUE + spectrum_output_file + bcolors.ENDC)
#
# ### Plotting the spectrum ###
# fig1 = plt.figure(figsize=(20, 10))
#
# font = {'family': 'serif',
#         'color': 'red',
#         'weight': 'normal',
#         'size': 15,
#         }
#
# # Setting plot labels
# plt.xlabel(r'$\lambda\,[\AA]$')
# plt.ylabel('Relative flux')
#
# # Remove sky spectrum from science spectrum
# spec = spec - skyspec
#
# # Setting plot ranges
# yoffset = (max(spec) - min(spec)) * 0.05
# pylab.ylim([min(spec) - yoffset, max(spec) + yoffset])
#
# plotoffset = (float(max(wavelengthrange)) - float(min(wavelengthrange))) * 0.01
# pylab.xlim(min(wavelengthrange) - plotoffset, max(wavelengthrange) + plotoffset)
#
# # Plot the actual data
# plt.plot(wavelengthrange, spec, 'b-')
# # plt.plot(wavelengthrange, spec - skyspec, 'b-')
# # plt.plot(wavelengthrange,spec2,'r-')
# # plt.show()
#
# # Setting plotpositions for ident lines
# if plot_idents and line_file != "":
#     plotminimum = min(spec) - yoffset
#     plotmaximum = max(spec) + yoffset
#     plotheigth = plotmaximum - plotminimum
#     plotmiddleplot = (plotminimum + plotmaximum) / 2.0
#     plotupperplot = plotminimum + 0.80 * plotheigth
#     plotlowerplot = plotminimum + 0.20 * plotheigth
#     plotuppercut1 = plotminimum + 0.70 * plotheigth
#     plotlowercut1 = plotminimum + 0.30 * plotheigth
#     plotuppercut2 = plotminimum + 0.68 * plotheigth
#     plotlowercut2 = plotminimum + 0.32 * plotheigth
#
#     # interpolate on data to find point for ident
#     f2 = interp1d(wavelengthrange, spec)
#     lines = open(line_file, "r")
#
#     # plot idents to figure
#     for line in lines:
#         liste = line.split()
#         if len(liste) == 1:
#             print(
#                 bcolors.WARNING + "     [WARNING] Broken identification found as '" + line + "', must consist of [wavelength(s) + name]. I will skip this one." + bcolors.ENDC)
#             continue
#         try:
#             float(liste[0])
#         except ValueError:
#             print(
#                 bcolors.WARNING + "     [WARNING] Broken identification found as '" + line + "', first entry not a number. I will skip this one." + bcolors.ENDC)
#             continue
#
#         #   Single ident plot
#         if len(liste) == 2:
#             ident_line = liste
#             wave_line = float(ident_line[0])
#             if float(wave_line) >= min(wavelengthrange) and float(wave_line) <= max(wavelengthrange):
#                 if f2(wave_line) <= plotmiddleplot:
#                     plt.plot(
#                         [wave_line, wave_line],
#                         [plotupperplot, f2(wave_line)],
#                         color='r',
#                         linestyle='-',
#                         linewidth=1.5,
#                     )
#                     plt.text(
#                         wave_line,
#                         plotupperplot,
#                         ident_line[1],
#                         rotation=90,
#                         ha='center',
#                         va='bottom',
#                         fontdict=font,
#                     )
#                 else:
#                     plt.plot(
#                         [wave_line, wave_line],
#                         [f2(wave_line), plotlowerplot],
#                         color='r',
#                         linestyle='-',
#                         linewidth=1.5,
#                     )
#                     plt.text(
#                         wave_line,
#                         plotlowerplot,
#                         ident_line[1],
#                         rotation=90,
#                         ha='center',
#                         va='top',
#                         fontdict=font,
#                     )
#         #   Multi ident plot
#         if len(liste) > 2:
#             points = []
#             for i in liste[:-1]:
#                 points.append(float(i))
#             pointcenter = sum(points) / float(len(points))
#             ident_name = str(liste[-1:])
#             # print(points," give ",pointcenter," bei ",liste[-1:])
#             if max(points) <= max(wavelengthrange) and min(points) >= min(wavelengthrange):
#                 if f2(pointcenter) <= plotmiddleplot:
#                     plt.plot(
#                         [pointcenter, pointcenter],
#                         [plotupperplot, plotuppercut1],
#                         color='r',
#                         linestyle='-',
#                         linewidth=1.5,
#                     )
#                     plt.text(
#                         pointcenter,
#                         plotupperplot,
#                         ident_name[2:-2],
#                         rotation=90,
#                         ha='center',
#                         va='bottom',
#                         fontdict=font,
#                     )
#                     for element in points:
#                         plt.plot(
#                             [element, element],
#                             [f2(element), plotuppercut2],
#                             color='r',
#                             linestyle='-',
#                             linewidth=1.5,
#                         )
#                         plt.plot(
#                             [pointcenter, element],
#                             [plotuppercut1, plotuppercut2],
#                             color='r',
#                             linestyle='-',
#                             linewidth=1.5,
#                         )
#                 if f2(pointcenter) > plotmiddleplot:
#                     plt.plot(
#                         [pointcenter, pointcenter],
#                         [plotlowerplot, plotlowercut1],
#                         color='r',
#                         linestyle='-',
#                         linewidth=1.5,
#                     )
#                     plt.text(
#                         pointcenter,
#                         plotlowerplot,
#                         ident_name[2:-2],
#                         rotation=90,
#                         ha='center',
#                         va='top',
#                         fontdict=font,
#                     )
#                     for element in points:
#                         plt.plot(
#                             [element, element],
#                             [f2(element), plotlowercut2],
#                             color='r',
#                             linestyle='-',
#                             linewidth=1.5,
#                         )
#                         plt.plot(
#                             [pointcenter, element],
#                             [plotlowercut1, plotlowercut2],
#                             color='r',
#                             linestyle='-',
#                             linewidth=1.5,
#                         )
#     lines.close()
#
# # Write the plot file
# plt.savefig(spectrum_output_file, bbox_inches='tight')
#
# plt.clf()
#
# ### Plotting the flatfield ###
# fig2 = plt.figure(figsize=(20, 10))
# # Setting plot labels
# plt.xlabel(r'$\lambda\,[\AA]$')
# plt.ylabel('Relative flux')
#
# # Setting plot ranges
# yoffset = (max(spec2) - min(spec2)) * 0.05
# pylab.ylim([min(spec2) - yoffset, max(spec2) + yoffset])
#
# plotoffset = (float(max(wavelengthrange)) - float(min(wavelengthrange))) * 0.01
# pylab.xlim(min(wavelengthrange) - plotoffset, max(wavelengthrange) + plotoffset)
#
# # Plot the actual data
# plt.plot(wavelengthrange, spec2, 'r-')
#
# print(bcolors.BOLD + "   Create flatfield plot " + bcolors.OKBLUE + 'flatfield.pdf' + bcolors.ENDC)
#
# # Write the plot file
# plt.savefig('flatfield.pdf', bbox_inches='tight')
#
# print(bcolors.BOLD + "   Write spectrum to file " + bcolors.OKBLUE + spectrum_output_data_file + bcolors.ENDC)
#
# #   Write spectrum table
# tbl = Table(names=['wave', 'flux', ], data=[wavelengthrange, spec, ])
# tbl.write(spectrum_output_data_file, format='ascii', overwrite=True)
#
# # os.system('touch '+spectrumdata)
# # specdatafile = open(spectrumdata,'w')
#
# # for i in range(len(wavelengthrange)):
# # specdatafile.write( str(wavelengthrange[i]) + "   " +str(spec[i]))
#
# print(bcolors.OKGREEN + "   Done" + bcolors.ENDC)
