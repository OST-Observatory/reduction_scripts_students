#! /home/pollux/.virtualenvs/photo/bin/python
# -*- coding: utf-8 -*-

'''
    First part of the reduction pipeline for data taken within the
    scope of the C7 observation of the astrophysics lab course at
    Potsdam University.

    All files can be given in one directory called 'rawfiles'. Alternatively
    images can be sorted into the following directory structure:
        * Images of the object
        * Dark frames
        * Flatfields
    If they are sorted into directories the FITS header keywords will be
    checked for consistency.

    Images in sub folders will be recognized, but only one level is considered.

   Version
   -------
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
'''

############################################################################
###            Configuration: modify the file in this section            ###
############################################################################

#######################  Simple folder structure  ##########################
rawfiles = '?'

##########################  Individual folders  ############################
### Path to the bias -- If set to '?', bias exposures are not used.
bias = '?'

### Path to the darks
darks = '?'

### Path to the flats
flats = '?'

### Path to the images
imgs  = '?'


############################################################################
###              Additional options: only edit if necessary              ###
############################################################################

#   Path to store the output (will usually be 'output',
#   but it can be changed as needed).
outdir = 'output/'

#   Remove cosmic rays?
rmcos = True
#rmcos = False

#   Tolerance between science and dark exposure times in s
tolerance = 4.


############################################################################
####                            Libraries                               ####
############################################################################

import tempfile

import warnings
warnings.filterwarnings('ignore')

from ost_photometry.reduce import redu
from ost_photometry.reduce import aux


############################################################################
###                                Main                                  ###
############################################################################

if __name__ == '__main__':
    ###
    #   Prepare directories and make checks
    #
    #   Create temporary directory
    temp_dir = tempfile.TemporaryDirectory()

    #   Prepare directories
    rawfiles = aux.prepare_reduction(
        outdir,
        bias,
        darks,
        flats,
        imgs,
        rawfiles,
        temp_dir,
        )


    ###
    #   Reduce images
    #
    redu.reduce_main(
        rawfiles,
        outdir,
        cosmics=rmcos,
        tolerance=tolerance,
        stack=False,
        shift_all=True,
        )
