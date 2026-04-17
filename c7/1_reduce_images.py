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
