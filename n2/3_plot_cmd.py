#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Small python script for CMD plotting

    TODO: replace the question marks ('?') with the adequate input

"""

############################################################################
#             Configuration: modify the file in this section               #
############################################################################

#   Name of the star cluster
name_of_star_cluster = "?"

###
#   Parameters regarding the file, containing the CMD data
#
#   Name of CMD data file
cmd_file_name = "?/cmd.dat"

#   Color and filter to plot in the CMD
#       -> the 'main_filter' defines the y-coordinate (ordinate) of the CMD
#       -> the second filter from the 'filter_list'
#           -> is used to calculate the color
#           -> defines the x-coordinate (abscissa) of the CMD
#           -> multiple entries in the list result in the creation of
#              multiple CMDs
main_filter = 'V'  # currently only 'V' is supported
filter_list = [
    # 'U',
    'B',
    # 'R',
    # 'I',
]

###
#   Calibration parameter
#
#   Calibration factor/zero point of the filter
cali = {
    'B': 0,
    'V': 0.,
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
x_plot_range_apparent = ["", ""]
y_plot_range_apparent = ["", ""]
#   Absolute CMD:
x_plot_range_absolute = ["", ""]
y_plot_range_absolute = ["", ""]

#   Size of the output figure in cm, default 8cm x 8cm
figure_size_x = "?"
figure_size_y = "?"

#   Name of output file, default: "cmd"
file_name = "cmd"

#   Filetype of output, supported filetypes:
#       -> png, pdf, ps, eps, and svg - default: pdf
file_type = "pdf"

#   Output directory
output_dir = "output/"

#   Plot error bars?
# do_error_bars = True
do_error_bars = False

############################################################################
#         Isochrone configuration: modify the file in this section         #
############################################################################

###
#   Isochrones in the archive
#   -> Set YAML file
#
#   NO isochrones
isochrone_configuration_file = ""

#   YY isochrones
isochrone_configuration_file = 'yy_isochrones.yaml'

#   basti-iac isochrones -> [Fe/H]=−1.58, Z = 0.0004, Y = 0.2476, [α/Fe]=0,
#   overshooting, diffusion, mass loss efficiency η = 0.3
isochrone_configuration_file = 'basti-iac_isochrones.yaml'

#   PARCES isochrones (CMD 3.6)
isochrone_configuration_file = 'parsec_3p6_isochrones.yaml'

#   PARCES isochrones (CMD 3.6, no TP-AGB evolution)
isochrone_configuration_file = 'parsec_3p6_noTP-AGB_isochrones.yaml'

############################################################################
#                               Libraries                                  #
############################################################################

import sys

import matplotlib

import numpy as np

from astropy.table import Table

from ost_photometry.analyze import plot, utilities
from ost_photometry.style import Bcolors

from ost_photometry import checks
from ost_photometry import utilities as base_utilities

############################################################################
#                           Routines & definitions                         #
############################################################################

matplotlib.rcParams['pdf.fonttype'] = 42

############################################################################
#                                  Main                                    #
############################################################################


if __name__ == '__main__':
    ###
    #   Check output directories
    #
    checks.check_output_directories(output_dir)

    ###
    #   Read CMD file
    #
    print(f'{Bcolors.BOLD}   Read file: {cmd_file_name}{Bcolors.ENDC}')

    #   Read table
    tbl_cmd = Table.read(cmd_file_name, format='ascii')
    if len(tbl_cmd) == 0:
        print(
            f'{Bcolors.FAIL}   The CMD table is empty => EXIT{Bcolors.ENDC}'
        )
        sys.exit()

    #   Check variables
    file_name, file_type = utilities.check_variable_apparent_cmd(
        file_name,
        file_type,
        main_filter,
        filter_list,
        cali,
    )

    #   Loop over all CMDs/colors
    for second_filer in filter_list:
        #   Extract data
        magnitude_filter_1 = tbl_cmd[f'{main_filter} [mag]'].value
        magnitude_filter_2 = tbl_cmd[f'{second_filer} [mag]'].value

        #   Apply zero point
        magnitude_filter_1 = magnitude_filter_1 + cali[main_filter]
        magnitude_filter_2 = magnitude_filter_2 + cali[second_filer]
        color = magnitude_filter_2 - magnitude_filter_1

        #   Get errors
        if do_error_bars:
            magnitude_filter_1_err = tbl_cmd[f'{main_filter}_err'].value
            magnitude_filter_2_err = tbl_cmd[f'{second_filer}_err'].value
            color_err = utilities.err_prop(
                magnitude_filter_1_err,
                magnitude_filter_2_err
            )
        else:
            magnitude_filter_1_err = None
            color_err = None

        ###
        #   Plot CMD
        #
        print(
            f'{Bcolors.BOLD}   Create {Bcolors.UNDERLINE} apparent'
            f'{Bcolors.ENDC}{Bcolors.BOLD}  CMD ({main_filter} vs. '
            f'{second_filer}-{main_filter}){Bcolors.ENDC}'
        )

        #   Plot apparent CMD
        plot.plot_apparent_cmd(
            color,
            magnitude_filter_1,
            name_of_star_cluster,
            file_name,
            file_type,
            main_filter,
            second_filer,
            figure_size_x=figure_size_x,
            figure_size_y=figure_size_y,
            y_plot_range_max=y_plot_range_apparent[1],
            y_plot_range_min=y_plot_range_apparent[0],
            x_plot_range_max=x_plot_range_apparent[1],
            x_plot_range_min=x_plot_range_apparent[0],
            output_dir=output_dir,
            magnitude_filter_1_err=magnitude_filter_1_err,
            color_err=color_err,
        )

        #   Check if the absolute CMD can be calculated
        if m_M == '?':
            if distance != '?':
                m_M = 5 * np.log10(float(distance) * 100.)
            else:
                m_M = 0.

        if eB_V != 0. and m_M != 0.:
            print('')
            print(
                f'\n{Bcolors.BOLD}   Create {Bcolors.UNDERLINE}absolute'
                f'{Bcolors.ENDC}{Bcolors.BOLD} CMD ({main_filter} vs. '
                f'{second_filer}-{main_filter}){Bcolors.ENDC}'
            )

            #   Read file with isochrone specification
            isochrone_configuration = base_utilities.read_params_from_yaml(
                isochrone_configuration_file
            )
            if isochrone_configuration:
                isochrones = isochrone_configuration.get('isochrones', '')
                isochrone_type = isochrone_configuration['isochrone_type']
                isochrone_column_type = isochrone_configuration['isochrone_column_type']
                isochrone_column = isochrone_configuration['isochrone_column']
                isochrone_keyword = isochrone_configuration['isochrone_keyword']
                isochrone_log_age = isochrone_configuration['isochrone_log_age']
                isochrone_legend = isochrone_configuration['isochrone_legend']

                #   Check isochrone parameters
                utilities.check_variable_absolute_cmd(
                    main_filter,
                    filter_list,
                    isochrone_column_type,
                    isochrone_column,
                )
            else:
                isochrones, isochrone_type, isochrone_column_type = '', '', ''
                isochrone_column, isochrone_keyword = '', ''
                isochrone_log_age, isochrone_legend = '', ''

            #   Correct for reddening and distance
            AV = RV * eB_V
            magnitude_filter_1 = magnitude_filter_1 - AV - m_M
            mag_color = color - eB_V

            #   Plot absolute CMD with isochrones
            plot.plot_absolute_cmd(
                mag_color,
                magnitude_filter_1,
                name_of_star_cluster,
                file_name,
                file_type,
                main_filter,
                second_filer,
                isochrones,
                isochrone_type,
                isochrone_column_type,
                isochrone_column,
                isochrone_log_age,
                isochrone_keyword,
                isochrone_legend,
                figure_size_x=figure_size_x,
                figure_size_y=figure_size_y,
                y_plot_range_max=y_plot_range_absolute[1],
                y_plot_range_min=y_plot_range_absolute[0],
                x_plot_range_max=x_plot_range_absolute[1],
                x_plot_range_min=x_plot_range_absolute[0],
                output_dir=output_dir,
                magnitude_filter_1_err=magnitude_filter_1_err,
                color_err=color_err,
            )

    print(f'{Bcolors.OKGREEN}   Done{Bcolors.ENDC}')

"""
    Change Log
    ----------
    * 20.11.2018
        - Change for automate plot range to trigger when non-float-able value
          given
        - Check separately if X and Y limitations are given or not and only use
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
"""
