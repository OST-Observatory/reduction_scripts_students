#! /home/pollux/.virtualenvs/photo/bin/python
# -*- coding: utf-8 -*-

############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

###
#   Name of the variable star
#
namestar     = "?"


###
#   Coordinates - Format:  ra = hh:mm:ss e.g. 19:44:42.8539591894
#                         dec = dd:am:as e.g. +54:49:42.887193554
#
ra_obj       = "??:??:??"
dec_obj      = "??:??:??"


###
#   Date of the minimum (UTC)
#   "yyyy:mm:ddThh:mm:ss" e.g., "2020-09-18T01:00:00"
#
transit_time = "yyyy:mm:ddThh:mm:ss"


###
#   Period (Algol: p=2.867315d, RZ Cas: p=1.1952499d, TV Cas: p=1.81259d)
#
period       = ?


####################  Program options  ######################

###
#   Finder options
#
#   Set sigma -> characterizes the size of the diffraction patterns
sigma = 3.0


############################################################################
####             Additional options: only edit if necessary             ####
####           Add up to 10 variable set for different filter           ####
############################################################################

########################  Filter 1  #########################
###
#   Define filter (e.g., U, B,V,...)
#
filter_1   = 'V'


###
#   Path to the images of filter 1
#
path_1     = './output/V/'


########################  Filter 2  #########################
###
#   Define filter (e.g., U, B,V,...)
#
filter_2   = 'B'


###
#   Path to the images of filter 2
#
path_2     = './output/B/'


########################  Filter 3  #########################
###
#   Define filter (e.g., U, B,V,...)
#
filter_3   = 'R'


###
#   Path to the images of filter 2
#
path_3     = './output/R/'


##############  Additional program options  #################

###
#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
#
outdir='output/'

###
#   Aperture or ePSF photometry
#
#   aper or PSF
photometry = 'APER'
#photometry = 'PSF'


###
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
#calib_method = 'UCAC4'


'''
The Full GSC2.3.2 Catalogue:
    Often not Johnson B & V magnitudes, but similar passbands. Only those
    with codes 3 and 4 are Johnson magnitudes.
'''
#calib_method = 'GSC2.3'


'''
    URAT1 Catalog (Zacharias+ 2015)
    B and V are also from the APASS survey.
'''
#calib_method = 'URAT1


'''
    NOMAD Catalog (Zacharias+ 2005):
    V: Photometric magnitude in Optical V band between 500 and 600 nm
    B: Photometric magnitude in Optical B band between 400 and 500 nm
'''
#calib_method = 'NOMAD'


'''
    Homogeneous Means in the UBV System (Mermilliod 1991):
    Johnson V-band & Johnson B-band magnitude (only B-V given)
'''
#calib_method = 'HMUBV'


'''
    Guide Star Photometric Catalog V2.4 (Bucciarelli+ 2001):
    Johnson B,V,R-band magnitude
'''
#calib_method = 'GSPC2.4'


'''
    AAVSO Photo. All Sky Surv. DR9(Henden+,2016):
    Johnson V-band & Johnson B-band magnitude
'''
calib_method = 'APASS'


'''
    Swift/UVOT Serendipitous Source Catalog (Yershov, 2015):
    AB magnitudes: U-AB, B-AB, V-AB
'''
#calib_method = 'Swift/UVOT'


'''
    XMM-OM Serendipitous Source Survey Catalogue (XMM-SUSS5.0)
    (Page+, 2021):
    AB magnitudes: UmAB, BmAB, VmAB
'''
#calib_method = 'XMM-OM'


'''
    Optical-UV-IR survey of North Celestial Cap (Gorbikov+, 2014):
    Johnson VRI magnitudes
'''
#calib_method = 'VRI-NCC'


'''
    The USNO-B1.0 Catalog (Monet+ 2003):
    BRI magnitudes
'''
#calib_method = 'USNO-B1.0'


#   Magnitude limit of the calibration stars
mag_range = (12., 15.)


###
#   Aperture options
#
#   Extraction radius stars in arcsec or pixel
rstars = 4.

#   Extraction radius background (inner and outer radii) in arcsec or pixel
rbg_in  = 5.
rbg_out = 7.

#   Unit
#r_unit = 'pixel'
r_unit = 'arcsec'


###
#   Correlation options
#
#   ID of the reference image
ref_ID = 0

#   Maximal separation between two objects in arcsec
seplimit = 5.

#   Limit for the number of images on which an object is not found.
#   When this limit is reached, the corresponding object is discarded.
nmissed = 10


###
#   Light curve options
#
#   Binning in days (set to None to deactivate)
tbin = 0.0001
tbin = 0.0002
tbin = None


############################################################################
####                            Libraries                               ####
############################################################################

import time

import warnings
warnings.filterwarnings('ignore')

from ost_photometry import style
from ost_photometry.analyze import analyze

import astropy.units as u

############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    #   Set start time
    start_time = time.time()

    #   Prepare variable lists and dictionaries from the individual
    #   definitions above
    filter_list = []
    img_paths = {}
    sigma_psf = {}
    for i in range(0,10):
        if 'filter_'+str(i) in locals():
            filter_list.append(locals()['filter_'+str(i)])
            img_paths[locals()['filter_'+str(i)]] = locals()['path_'+str(i)]
            sigma_psf[locals()['filter_'+str(i)]] = sigma


    ###
    #   Initialize image ensemble container
    #
    img_container = analyze.image_container()


    ###
    #   Extract flux
    #
    analyze.extract_flux_multi(
        img_container,
        filter_list,
        namestar,
        img_paths,
        outdir,
        sigma_psf,
        ra_obj,
        dec_obj,
        photometry=photometry,
        rstars=rstars,
        rbg_in=rbg_in,
        rbg_out=rbg_out,
        r_unit=r_unit,
        ref_ID=ref_ID,
        nmissed=nmissed,
        seplimit=seplimit*u.arcsec,
        #correl_method='own',
        )


    ###
    #   Calibrate data and plot light curves
    analyze.calibrate_data_mk_lc(
        img_container,
        filter_list,
        ra_obj,
        dec_obj,
        namestar,
        outdir,
        transit_time,
        period,
        binn=tbin,
        ref_ID=ref_ID,
        calib_method=calib_method,
        nmissed=nmissed,
        photo_type=photometry,
        seplimit=seplimit*u.arcsec,
        #correl_method='own',
        plot_sigma=True,
        #derive_Tcs=True,
        )


print(style.bcolors.OKGREEN+"   Done"+style.bcolors.ENDC)
print("--- %s minutes ---" % ((time.time() - start_time)/60.))
