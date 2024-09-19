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

#   List of colors to be plotted in the CMDs. The color filters must be
#   specified with a dash between them, such as B-V. The second filter will be
#   the magnitude plotted on the ordinate while the color is plotted on the
#   abscissa. (B-V) means that V is plotted against B-V. Multiple entries in
#   the list result in the generation of multiple CMDs.
filter_color_combinations = [
    "B-V",
]

###
#   Calibration parameter
#
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

#   X and Y range to plot (change according to your data)
#       - Plot range is automatically adjusted if range is set to ""
#       - Example for 2 CMD plots aka 3 available filters:
#           x_plot_range_apparent = [(0., 2.), (-1., 1.)]
#   Apparent CMD:
x_plot_range_apparent = [("", ""), ]
y_plot_range_apparent = [("", ""), ]
#   Absolute CMD:
x_plot_range_absolute = [("", ""), ]
y_plot_range_absolute = [("", ""), ]

#   Size of the output figure in cm, default 8cm x 8cm
figure_size_x = "?"
figure_size_y = "?"

#   Name of output file, default: "cmd"
file_name = "cmd"

#   Filetype of output, supported filetypes:
#       -> png, pdf, ps, eps, and svg - default: pdf
file_type = "pdf"

#   Output directory
output_dir = "output"

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

from ost_photometry.analyze import plots, utilities
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

    #   Loop over all CMDs/colors
    for filter_id, color in enumerate(filter_color_combinations):
        filter_list = color.split('-')
        filter_1 = filter_list[0]
        filter_2 = filter_list[1]

        #   Check variables
        file_name, file_type = utilities.check_variable_apparent_cmd(
            file_name,
            file_type,
        )

        #   Extract data
        magnitude_filter_1 = tbl_cmd[f'{filter_1} [mag]'].value
        magnitude_filter_2 = tbl_cmd[f'{filter_2} [mag]'].value

        #   Calculate color
        color = magnitude_filter_1 - magnitude_filter_2

        #   Get errors
        if do_error_bars:
            magnitude_filter_1_err = tbl_cmd[f'{filter_1}_err'].value
            magnitude_filter_2_err = tbl_cmd[f'{filter_2}_err'].value
            color_err = utilities.err_prop(
                magnitude_filter_1_err,
                magnitude_filter_2_err
            )
        else:
            magnitude_filter_2_err = None
            color_err = None

        ###
        #   Plot CMD
        #
        print(
            f'{Bcolors.BOLD}   Create {Bcolors.UNDERLINE}apparent'
            f'{Bcolors.ENDC}{Bcolors.BOLD} CMD: {filter_2} vs. '
            f'{filter_1}-{filter_2}{Bcolors.ENDC}'
        )

        #   Setup CMD object
        cmds = plots.MakeCMDs(
            name_of_star_cluster,
            file_name,
            file_type,
            filter_2,
            filter_1,
            color,
            magnitude_filter_2,
            color_err=color_err,
            magnitude_filter_2_err=magnitude_filter_2_err,
            output_dir=output_dir,
        )

        #   Plot apparent CMD
        cmds.plot_apparent_cmd(
            figure_size_x=figure_size_x,
            figure_size_y=figure_size_y,
            y_plot_range_max=y_plot_range_apparent[filter_id][1],
            y_plot_range_min=y_plot_range_apparent[filter_id][0],
            x_plot_range_max=x_plot_range_apparent[filter_id][1],
            x_plot_range_min=x_plot_range_apparent[filter_id][0],
        )

        #   Check if the absolute CMD can be calculated
        if m_M == '?':
            if distance != '?':
                m_M = 5 * np.log10(float(distance) * 100.)
            else:
                m_M = 0.

        if m_M != 0.:
            print(
                f'{Bcolors.BOLD}   Create {Bcolors.UNDERLINE}absolute'
                f'{Bcolors.ENDC}{Bcolors.BOLD} CMD: {filter_2} vs. '
                f'{filter_1}-{filter_2}{Bcolors.ENDC}'
            )

            #   Read file with isochrone specification
            isochrone_configuration = base_utilities.read_params_from_yaml(isochrone_configuration_file)
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
                    filter_list,
                    isochrone_column_type,
                    isochrone_column,
                )
            else:
                isochrones, isochrone_type, isochrone_column_type = '', '', ''
                isochrone_column, isochrone_keyword = '', ''
                isochrone_log_age, isochrone_legend = '', ''

            #   Plot absolute CMD with isochrones
            cmds.plot_absolute_cmd(
                eB_V,
                m_M,
                isochrones,
                isochrone_type,
                isochrone_column_type,
                isochrone_column,
                isochrone_log_age,
                isochrone_keyword,
                isochrone_legend,
                rv=RV,
                figure_size_x=figure_size_x,
                figure_size_y=figure_size_y,
                y_plot_range_max=y_plot_range_absolute[filter_id][1],
                y_plot_range_min=y_plot_range_absolute[filter_id][0],
                x_plot_range_max=x_plot_range_absolute[filter_id][1],
                x_plot_range_min=x_plot_range_absolute[filter_id][0],
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
    * 13.03.2024
        - Switch from function-based CMD plotting to class-based plotting
        - Removed calibration parameter
"""
