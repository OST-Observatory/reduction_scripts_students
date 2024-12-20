#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    First part of the reduction pipeline for data taken within the
    scope of the C7 observation of the astrophysics lab course at
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
#                          Simple folder structure                         #
############################################################################
raw_files: str = '?'

############################################################################
#                              Individual folders                          #
############################################################################
# Path to the bias -- If set to '?', bias exposures are not used.
bias: str = '?'

# Path to the darks
darks: str = '?'

# Path to the flats
flats: str = '?'

# Path to the images
images: str = '?'


############################################################################
#                Additional options: only edit if necessary                #
############################################################################

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir: str = 'output/'

#   Remove cosmic rays?
rm_cosmic_rays: bool = True
# rm_cosmic_rays: bool = False

#   Tolerance between science and dark exposure times in s
exposure_time_tolerance: float = 5.

#   Tolerance between the camera chip temperatures of the images
temperature_tolerance: float = 5.

#   Number of cores used for multiprocessing
n_cores_multiprocessing: int = 4

############################################################################
#                               Libraries                                  #
############################################################################

import time

import tempfile

import warnings
warnings.filterwarnings('ignore')

from astropy import log
log.setLevel('ERROR')

from ost_photometry.reduce import redu
from ost_photometry.reduce import utilities
from ost_photometry import style


############################################################################
#                                  Main                                    #
############################################################################

if __name__ == '__main__':
    #   Set start time
    start_time = time.time()

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
        stack_images=False,
        shift_all=True,
        temperature_tolerance=temperature_tolerance,
        n_cores_multiprocessing=n_cores_multiprocessing,
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.))

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