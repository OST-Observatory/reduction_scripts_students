#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    First part of the reduction pipeline for data taken within the
    scope of the N2 observation of the astrophysics lab course at
    Potsdam University.

    All files can be given in one directory called 'raw_files'. Alternatively
    images can be sorted into the following directory structure:
        * Images of the object
        * Dark frames
        * Flatfields
    If they are sorted into directories the FITS header keywords will be
    checked for consistency.

    Images in sub folders will be recognized, but only one level is considered.
"""

############################################################################
#                         Simple folder structure                          #
############################################################################
raw_files = '?'

############################################################################
#                           Individual folders                             #
############################################################################
# Path to the bias -- If set to '?', bias exposures are not used.
bias = '?'

# Path to the darks
darks = '?'

# Path to the flats
flats = '?'

# Path to the images
images = '?'


############################################################################
#                Additional options: only edit if necessary                #
############################################################################

#   Cluster identifier (e.g., NGC 381):
#   Only set target if you want to filter the cluster images by the target
#   name. For this to work, the target name must appear as a FITS header
#   keyword.
target_name = None

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir = 'output/'

#   Remove cosmic rays?
rm_cosmic_rays = True
# rm_cosmic_rays = False

#   Tolerance between science and dark exposure times in s
exposure_time_tolerance = 4.


############################################################################
#                               Libraries                                  #
############################################################################

import tempfile

import warnings
warnings.filterwarnings('ignore')

from astropy import log
log.setLevel('ERROR')

from ost_photometry.reduce import redu
from ost_photometry.reduce import utilities


############################################################################
#                                  Main                                    #
############################################################################

if __name__ == '__main__':
    ###
    #   Prepare directories and make checks
    #
    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Prepare directories
    raw_files = utilities.prepare_reduction(
        output_dir,
        bias,
        darks,
        flats,
        images,
        raw_files,
        temp_dir,
        )

    ###
    #   Reduce images
    #
    redu.reduce_main(
        raw_files,
        output_dir,
        rm_cosmic_rays=rm_cosmic_rays,
        exposure_time_tolerance=exposure_time_tolerance,
        target_name=target_name,
        )

"""
    Change Log
    ----------
        0.1   (18.11.2020)
           - initial release
        0.2   (12.01.2021)
           - almost complete rewrite
           - switched to lacosmics from ccdproc for cosmics removal
           - added flexibility
           - added support for an arbitrary number of filters
           - reduced number of code lines
        0.21  (13.01.2021)
           - fixed header handling
           - further flexibility updates
        0.22 (28.01.2021)
           - minor style updates
        0.23 (30.06.2021)
           - minor style update
        0.3  (10.02.2022)
           - complete rewrite using ccdproc
        0.4  (07.09.2022)
           - set a couple of default parameters
        0.5  (28.08.2023)
           - some adjustments because of pipeline updates
"""