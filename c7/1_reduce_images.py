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

#   Number of cores used for multiprocessing (``None`` = half the CPUs)
n_cores_multiprocessing: int | None = None

#   Measure the quality of every reduced frame (FWHM, roundness, number of
#   stars, sky background). Results: ``output/frame_quality.ecsv``, the FITS
#   headers, and ``output/diagnostics/frame_quality/``. Needed for the two
#   options below.
measure_frame_quality: bool = True

#   Reject poor frames before the alignment. ``None`` keeps every frame
#   (default for light curves: a rejected frame is a missing epoch). Example
#   that drops frames with a strongly deviating FWHM or almost no stars
#   (clouds); rejected frames are moved to ``output/rejected_lights/``:
frame_selection: dict | None = None
# frame_selection: dict | None = {"fwhm_sigma_clip": 3.0, "n_stars_min": 20}

#   Alignment reference: ``best_fwhm`` = sharpest frame (per filter), or
#   ``first`` = first frame in time (the old behaviour).
reference_image_selection: str = "best_fwhm"

#   Alignment backend. ``aa_true`` (astroalign similarity) is the C7 default.
#   ``wcs`` reprojects onto a shared sky grid (handles rotation / flip);
#   needs a current ost_photometry so the bad-pixel mask is not smeared.
shift_method: str = "aa_true"
# shift_method: str = "wcs"

#   Warp every science frame onto the reference grid. Set False for native
#   pixels; then ``2_obtain_flux.py`` must keep ``wcs_solve_all_images=True``
#   so analysis matches on the sky.
shift_all: bool = True
# shift_all: bool = False

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
        shift_all=shift_all,
        shift_method=shift_method,
        temperature_tolerance=temperature_tolerance,
        n_cores_multiprocessing=n_cores_multiprocessing,
        measure_frame_quality=measure_frame_quality,
        frame_selection=frame_selection,
        reference_image_selection=reference_image_selection,
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.))
