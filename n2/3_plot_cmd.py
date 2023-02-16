#! /usr/bin/python
# -*- coding: utf-8 -*-

'''
    Small python script for CMD plotting

    TODO: replace the question marks ('?') with the adequate input

#### Change Log

* 20.11.2018
    - Change for automate plot range to trigger when non-floatable value
      given
    - Check separately if X and Y limitaions are given or not and only use
      for automated where non is given
    - Switched to TrueType fonts for the pdf
    - Added optional labels for Isochrones
* 20.11.2020
    - Added comments
    - Added line style cycler
    - Added option to calibrate the CMD
* 18.04.2021
    - Almost complete rewrite
    - Added support for simultaneous calculation of apparent and absolute
      CMD
    - Added support for isochrone files that contain several isochrones
* 31.08.2022
    - Add error bars
    - Add additional isochrones
'''

############################################################################
####          Configuration: modify the file in this section            ####
############################################################################

#   Name of the star cluster
nameOfStarcluster = "?"

###
#   Parameters regarding the file, containing the CMD data
#
#   Name of CMD data file
CMDFileName = "?/cmd.dat"

#   Color and filter to plot in the CMD
#       -> 'filt_1' defines the y-coordinate (ordinate) of the CMD
#       -> 'filt_2' second filter used to calculate the color
#           -> Defines the x-coordinate (abscissa) of the CMD
#           -> Multiple entries result in the creation of multiple CMDs
filt_1 = 'V'                     #  currently only 'V' is supported
filt_2 = [
    #'U',
    'B',
    #'R',
    #'I',
    ]


###
#   Calibration parameter
#
#   Calibration factor/zero point of the filter
cali = {
    'B':0,
    'V':0.,
    }

#   E_B-V of the cluster
eB_V = 0.

#   R_V
RV = 3.1

#   Give either distance modulus of the cluster or the distance in kpc
m_M = '?'

distance = '?'



###
#   Plot parameter
#

#   x_Range=[xRangeMin:xRangeMax] & y_Range=[yRangeMin:yRangeMax]
#       -> x and y range to plot (change according to your data)
#       -> The plot range is automatically adjusted, if range is set to ""
#   Apparent CMD:
x_Range_apparent=["",""]
y_Range_apparent=["",""]
#   Absolute CMD:
x_Range_absolute=["",""]
y_Range_absolute=["",""]

#   Size of the output figure in cm, default 8cm x 8cm
size_x = "?"
size_y = "?"

#   Name of output file, default: "cmd"
filename    = "cmd"

#   Filetype of output, supported filetypes:
#       -> png, pdf, ps, eps, and svg - default: pdf
filetype    = "pdf"

#   Output directory
outdir = "output/"

#   Plot error bars?
#do_error_bars = True
do_error_bars = False


############################################################################
####      Isochrone configuration: modify the file in this section      ####
############################################################################

###
#   Isochrones in the archive
#   -> Set YAML file
#
#   NO isochrones
iso_config_file = ""

#   YY isochrones
iso_config_file = 'yy_isochrones.yaml'

#   basti-iac isochrones -> [Fe/H]=−1.58, Z = 0.0004, Y = 0.2476, [α/Fe]=0,
#   overshooting, diffusion, mass loss efficiency η = 0.3
iso_config_file = 'basti-iac_isochrones.yaml'

#   PARCES isochrones (CMD 3.6)
iso_config_file = 'parsec_3p6_isochrones.yaml'

#   PARCES isochrones (CMD 3.6, no TP-AGB evolution)
iso_config_file = 'parsec_3p6_noTP-AGB_isochrones.yaml'


############################################################################
####                            Libraries                               ####
############################################################################

import sys

import matplotlib.pyplot as plt
import pylab, os, matplotlib

import numpy as np

from astropy.table import Table

from ost_photometry.analyze import plot, aux
from ost_photometry.style import bcolors

from ost_photometry import checks
from ost_photometry import aux as aux_general

############################################################################
####                        Routines & definitions                      ####
############################################################################

matplotlib.rcParams['pdf.fonttype'] = 42


############################################################################
####                               Main                                 ####
############################################################################


if __name__ == '__main__':
    ###
    #   Check output directories
    #
    checks.check_out(outdir)


    ###
    #   Read CMD file
    #
    print(f'{bcolors.BOLD}   Read file: {CMDFileName}{bcolors.ENDC}')

    #   Read table
    tbl_cmd = Table.read(CMDFileName, format='ascii')
    if len(tbl_cmd) == 0:
        print(
            f'{bcolors.FAIL}   The CMD table is empty => EXIT{bcolors.ENDC}'
            )
        sys.exit()

    #   Read file with isochrone specification
    iso_config = aux_general.read_params_from_yaml(iso_config_file)
    isos = iso_config.get('isos', '')
    isotype = iso_config['isotype']
    ISOcolumntype = iso_config['ISOcolumntype']
    ISOcolumn = iso_config['ISOcolumn']
    keyword = iso_config['keyword']
    logAGE = iso_config['logAGE']
    IsoLabels = iso_config['IsoLabels']

    #   Check variable
    filename, filetype = aux.check_variable(
        filename,
        filetype,
        filt_1,
        filt_2,
        cali,
        ISOcolumntype,
        ISOcolumn,
        )

    #   Loop over all CMDs/colors
    for fil in filt_2:
        #   Set color
        color = fil+'-'+filt_1

        #   Extract data
        #mag_filt_1 = tbl_cmd[filt_1+' [mag] (0)'].value
        #mag_filt_2 = tbl_cmd[fil+' [mag] (0)'].value
        mag_filt_1 = tbl_cmd[filt_1+' [mag]'].value
        mag_filt_2 = tbl_cmd[fil+' [mag]'].value

        #   Apply zero point
        mag_filt_1  = mag_filt_1 + cali[filt_1]
        mag_filt_2  = mag_filt_2 + cali[fil]
        color = mag_filt_2 - mag_filt_1

        #   Get errors
        if do_error_bars:
            #mag_filt_1_err = tbl_cmd[filt_1+'_err [mag] (0)'].value
            #mag_filt_2_err = tbl_cmd[fil+'_err [mag] (0)'].value
            mag_filt_1_err = tbl_cmd[filt_1+'_err [mag]'].value
            mag_filt_2_err = tbl_cmd[fil+'_err [mag]'].value
            color_err = aux.err_prop(mag_filt_1_err, mag_filt_2_err)
        else:
            mag_filt_1_err = None
            color_err = None

        ###
        #   Plot CMD
        #
        print(
            f'{bcolors.BOLD}   Create {bcolors.UNDERLINE} apparent'
            f'{bcolors.ENDC}{bcolors.BOLD}  CMD ({filt_1} vs. {fil}-'
            f'{filt_1}){bcolors.ENDC}'
            )

        #   Plot apparent CMD
        plot.plot_apparent_cmd(
            color,
            mag_filt_1,
            nameOfStarcluster,
            filename,
            filetype,
            filt_1,
            fil,
            size_x=size_x,
            size_y=size_y,
            yRangeMax=y_Range_apparent[1],
            yRangeMin=y_Range_apparent[0],
            xRangeMax=x_Range_apparent[1],
            xRangeMin=x_Range_apparent[0],
            outdir=outdir,
            mag_filt_err=mag_filt_1_err,
            color_err=color_err,
            )

        #   Check if the absolute CMD can be calculated
        if m_M == '?':
            if distance != '?':
                m_M = 5 * np.log10(float(distance)*100.)
            else:
                m_M = 0.

        if eB_V != 0. and m_M != 0.:
            print('')
            print(
                f'\n{bcolors.BOLD}   Create {bcolors.UNDERLINE}absolute'
                f'{bcolors.ENDC}{bcolors.BOLD} CMD ({filt_1} vs. {fil}-'
                f'{filt_1}){bcolors.ENDC}'
                )

            #   Correct for reddening and distance
            AV        = RV*eB_V
            mag_filt_1  = mag_filt_1-AV-m_M
            mag_color = color-eB_V

            #   Plot absolute CMD with isochrones
            plot.plot_absolute_cmd(
                mag_color,
                mag_filt_1,
                nameOfStarcluster,
                filename,
                filetype,
                filt_1,
                fil,
                isos,
                isotype,
                ISOcolumntype,
                ISOcolumn,
                logAGE,
                keyword,
                IsoLabels,
                size_x=size_x,
                size_y=size_x,
                yRangeMax=y_Range_absolute[1],
                yRangeMin=y_Range_absolute[0],
                xRangeMax=x_Range_absolute[1],
                xRangeMin=x_Range_absolute[0],
                outdir=outdir,
                mag_filt_err=mag_filt_1_err,
                color_err=color_err,
                )

    print(f'{bcolors.OKGREEN}   Done{bcolors.ENDC}')

