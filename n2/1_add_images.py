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
#   PLEASE NOTE: Either specify a single directory with all raw files at
#                this point, or use the directory structure shown below.
raw_files: str = "?"

############################################################################
#                           Individual folders                             #
############################################################################
#   PLEASE NOTE: Either use the directory scheme shown below or specify a
#   single directory with all raw files above.

# Path to the bias -- If set to '?', bias exposures are not used.
bias: str = "?"

# Path to the darks
darks: str = "?"

# Path to the flats
flats: str = "?"

# Path to the images
images: str = "?"


############################################################################
#                Additional options: only edit if necessary                #
############################################################################

#   Cluster identifier (e.g., NGC 381):
#   Only set target if you want to filter the cluster images by the target
#   name. For this to work, the target name must appear as a FITS header
#   keyword.
target_name: str | None = None

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
output_dir: str = "output/"

#   Remove cosmic rays?
rm_cosmic_rays: bool = True
# rm_cosmic_rays: bool = False

#   Tolerance between science and dark exposure times in s
exposure_time_tolerance: float = 5.0

#   Tolerance between the camera chip temperatures of the images
temperature_tolerance: float = 5.0

#   Number of cores used for multiprocessing
n_cores_multiprocessing: int = 8


############################################################################
#                               Libraries                                  #
############################################################################

import tempfile
import time
import warnings

warnings.filterwarnings("ignore")

from astropy import log

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.reduce import redu, utilities

############################################################################
#                                  Main                                    #
############################################################################

if __name__ == "__main__":
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
        target_name=target_name,
        temperature_tolerance=temperature_tolerance,
        n_cores_multiprocessing=n_cores_multiprocessing,
        # Pixel match (default) or WCS reproject: shift_method="wcs"
        # (uses wcs_method, default ASTAP, if a frame has no celestial WCS).
        shift_method="aa_true",
        shift_all=True,
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
