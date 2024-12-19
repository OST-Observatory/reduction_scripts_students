#! /usr/bin/python
# -*- coding: utf-8 -*-

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
# calibration_method: str = 'UCAC4'


'''
The Full GSC2.3.2 Catalogue:
    Often not Johnson B & V magnitudes, but similar passbands. Only those
    with codes 3 and 4 are Johnson magnitudes.
'''
# calibration_method: str = 'GSC2.3'


'''
    URAT1 Catalog (Zacharias+ 2015)
    B and V are also from the APASS survey.
'''
# calibration_method: str = 'URAT1


'''
    NOMAD Catalog (Zacharias+ 2005):
    V: Photometric magnitude in Optical V band between 500 and 600 nm
    B: Photometric magnitude in Optical B band between 400 and 500 nm
'''
# calibration_method: str = 'NOMAD'


'''
    Homogeneous Means in the UBV System (Mermilliod 1991):
    Johnson V-band & Johnson B-band magnitude (only B-V given)
'''
# calibration_method: str = 'HMUBV'


'''
    Guide Star Photometric Catalog V2.4 (Bucciarelli+ 2001):
    Johnson B,V,R-band magnitude
'''
# calibration_method: str = 'GSPC2.4'


'''
    AAVSO Photo. All Sky Surv. DR9(Henden+,2016):
    Johnson V-band & Johnson B-band magnitude
'''
calibration_method: str = 'APASS'

'''
    Swift/UVOT Serendipitous Source Catalog (Yershov, 2015):
    AB magnitudes: U-AB, B-AB, V-AB
'''
# calibration_method: str = 'Swift/UVOT'


'''
    XMM-OM Serendipitous Source Survey Catalogue (XMM-SUSS5.0)
    (Page+, 2021):
    AB magnitudes: UmAB, BmAB, VmAB
'''
# calibration_method: str = 'XMM-OM'


'''
    Optical-UV-IR survey of North Celestial Cap (Gorbikov+, 2014):
    Johnson VRI magnitudes
'''
# calibration_method: str = 'VRI-NCC'


'''
    The USNO-B1.0 Catalog (Monet+ 2003):
    BRI magnitudes
'''
# calibration_method: str = 'USNO-B1.0'


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
reference_image_id: int = 0

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

    ###
    #   Initialize observation container
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

    ###
    #   Extract flux
    #
    observation.extract_flux_multi(
        filter_list,
        image_paths,
        output_dir,
        fwhm_object_psf=fwhm_object_psf,
        photometry_extraction_method=photometry_extraction_method,
        radius_aperture=radius_aperture,
        inner_annulus_radius=inner_annulus_radius,
        outer_annulus_radius=outer_annulus_radius,
        radii_unit=radii_unit,
        reference_image_id=reference_image_id,
        n_allowed_non_detections_object=n_allowed_non_detections_object,
        separation_limit=separation_limit * u.arcsec,
        # duplicate_handling_object_identification={'astropy': 'flux'},
        # correlation_method='own',
    )

    ###
    #   Calibrate data and plot light curves
    observation.calibrate_data_mk_light_curve(
        filter_list,
        output_dir,
        binning_factor=binning_factor,
        calibration_method=calibration_method,
        n_allowed_non_detections_object=n_allowed_non_detections_object,
        photometry_extraction_method=photometry_extraction_method,
        separation_limit=separation_limit * u.arcsec,
        # correlation_method='own',
        plot_sigma=True,
        # derive_transformation_coefficients=True,
        calculate_zero_point_statistic=False,
        # duplicate_handling_object_identification={'astropy': 'flux'},
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.))
