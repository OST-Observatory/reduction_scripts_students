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

#   Remove cosmic rays (L.A.Cosmic)? ``"auto"`` leaves them to the stack for
#   filters with at least ``cosmic_ray_auto_min_frames`` images (the
#   sigma-clipped stack rejects them there) and removes them otherwise.
rm_cosmic_rays: bool | str = "auto"
# rm_cosmic_rays: bool | str = True
# rm_cosmic_rays: bool | str = False
cosmic_ray_auto_min_frames: int = 7

#   Gain and read noise of the camera: ``"auto"`` (measured from the bias and
#   flat frames when the FITS header does not name the readout mode, as with
#   MaxIm DL's "Fast" / "Normal", else catalog), ``"catalog"`` (FITS header /
#   camera catalog) or ``"measured"`` (always from the bias and flat frames
#   of this data set; falls back to the catalog if that is not possible).
camera_noise_source: str = "auto"
# camera_noise_source: str = "measured"

#   Tolerance between science and dark exposure times in s
exposure_time_tolerance: float = 5.0

#   Tolerance between the camera chip temperatures of the images
temperature_tolerance: float = 5.0

#   Number of cores used for multiprocessing (``None`` = half the CPUs)
n_cores_multiprocessing: int | None = None

#   Alignment backend. ``aa_true`` (astroalign similarity; pixel match) is
#   the default. ``wcs`` reprojects onto a shared sky grid (uses
#   ``wcs_method``, default ASTAP, if a frame has no celestial WCS).
shift_method: str = "aa_true"
# shift_method: str = "wcs"

#   Warp every science frame onto one common reference grid (all filters).
shift_all: bool = True

#   Measure the quality of every reduced frame (FWHM, roundness, number of
#   stars, sky background). Results: ``output/frame_quality.ecsv``, the FITS
#   headers, and ``output/diagnostics/frame_quality/``. Needed for the
#   selection, the weighting and the ``best_fwhm`` reference below.
measure_frame_quality: bool = True

#   Reject poor frames before stacking. ``None`` stacks every frame. Example:
#   keep the sharpest 90 % per filter (frames without stars are dropped as
#   well); rejected frames are moved to ``output/rejected_lights/``:
frame_selection: dict | None = None
# frame_selection: dict | None = {"best_fraction": 0.9}

#   Per-frame weights in the stack: ``none`` (equal weights), ``fwhm``,
#   ``n_stars`` or ``noise``. Weighting sharpens the stack but changes its
#   effective PSF and noise; for photometry on the stack ``none`` is the
#   conservative choice.
stack_weighting: str = "none"

#   Alignment reference: ``best_fwhm`` = sharpest frame of the data set, or
#   ``first`` = first frame in time (the old behaviour).
reference_image_selection: str = "best_fwhm"

#   Keep the registered single frames (``output/aligned_lights/``) after
#   stacking, e.g. to re-stack with another selection.
keep_aligned_lights: bool = False


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
        cosmic_ray_auto_min_frames=cosmic_ray_auto_min_frames,
        camera_noise_source=camera_noise_source,
        exposure_time_tolerance=exposure_time_tolerance,
        target_name=target_name,
        temperature_tolerance=temperature_tolerance,
        n_cores_multiprocessing=n_cores_multiprocessing,
        shift_method=shift_method,
        shift_all=shift_all,
        measure_frame_quality=measure_frame_quality,
        frame_selection=frame_selection,
        stack_weighting=stack_weighting,
        reference_image_selection=reference_image_selection,
        keep_aligned_lights=keep_aligned_lights,
    )

    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
