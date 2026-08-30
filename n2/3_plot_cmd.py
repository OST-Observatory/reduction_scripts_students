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
name_of_star_cluster: str = "?"

###
#   Parameters regarding the file, containing the CMD data
#
#   Pipeline post-processed table (epoch-native ECSV) or legacy ASCII

#   Name of CMD data file
cmd_file_name: str = "?/cmd.dat"

#   List of colors to be plotted in the CMDs. The color filters must be
#   specified with a dash between them, such as B-V. The second filter will be
#   the magnitude plotted on the ordinate while the color is plotted on the
#   abscissa. (B-V) means that V is plotted against B-V. Multiple entries in
#   the list result in the generation of multiple CMDs.
filter_color_combinations: list[str] = [
    "B-V",
]

###
#   Calibration parameter
#
#   E_B-V of the cluster
eB_V: float = 0.

#   R_V
RV: float = 3.1

#   Give either distance modulus of the cluster or the distance in kpc
m_M: str | float = '?'

distance: str | float = '?'

###
#   Plot parameter
#

#   X and Y range to plot (change according to your data)
#       - Plot range is automatically adjusted if range is set to ""
#       - Example for 2 CMD plots aka 3 available filters:
#           x_plot_range_apparent = [(0., 2.), (-1., 1.)]
#   Apparent CMD:
x_plot_range_apparent: list[tuple[str | float, str | float]] = [("", ""), ]
y_plot_range_apparent: list[tuple[str | float, str | float]] = [("", ""), ]
#   Absolute CMD:
x_plot_range_absolute: list[tuple[str | float, str | float]] = [("", ""), ]
y_plot_range_absolute: list[tuple[str | float, str | float]] = [("", ""), ]

#   Size of the output figure in cm, default 8cm x 8cm
figure_size_x: str | float = "?"
figure_size_y: str | float = "?"

#   Name of output file, default: "cmd"
file_name: str = "cmd"

#   Filetype of output, supported filetypes:
#       -> png, pdf, ps, eps, and svg - default: pdf
file_type: str = "pdf"

#   Output directory
output_dir: str = "output"

#   Plot error bars?
# do_error_bars: bool = True
do_error_bars: bool = False

############################################################################
#         Isochrone configuration: modify the file in this section         #
############################################################################

###
#   Isochrones in the archive
#   -> Set YAML file
#
#   NO isochrones
isochrone_configuration_file: str = ""

#   YY isochrones
isochrone_configuration_file: str = 'yy_isochrones.yaml'

#   basti-iac isochrones -> [Fe/H]=−1.58, Z = 0.0004, Y = 0.2476, [α/Fe]=0,
#   overshooting, diffusion, mass loss efficiency η = 0.3
isochrone_configuration_file: str = 'basti-iac_isochrones.yaml'

#   PARCES isochrones (CMD 3.6)
isochrone_configuration_file: str = 'parsec_3p6_isochrones.yaml'

#   PARCES isochrones (CMD 3.6, no TP-AGB evolution)
isochrone_configuration_file: str = 'parsec_3p6_noTP-AGB_isochrones.yaml'

############################################################################
#                               Libraries                                  #
############################################################################

import sys

import matplotlib
from ost_photometry import checks
from ost_photometry.analyze.cmd_prepare import load_cmd_table
from ost_photometry.analyze.plots import plot_cmds_from_table
from ost_photometry.style import Bcolors

matplotlib.rcParams['pdf.fonttype']: int = 42

############################################################################
#                                  Main                                    #
############################################################################


if __name__ == '__main__':
    checks.check_output_directories(output_dir)

    print(f'{Bcolors.BOLD}   Read file: {cmd_file_name}{Bcolors.ENDC}')
    tbl_cmd = load_cmd_table(cmd_file_name)
    if len(tbl_cmd) == 0:
        print(
            f'{Bcolors.FAIL}   The CMD table is empty => EXIT{Bcolors.ENDC}'
        )
        sys.exit()

    plot_cmds_from_table(
        tbl_cmd,
        filter_color_combinations,
        name_of_star_cluster=name_of_star_cluster,
        file_name=file_name,
        file_type=file_type,
        output_dir=output_dir,
        e_b_v=eB_V,
        rv=RV,
        m_m=m_M,
        distance=distance,
        do_error_bars=do_error_bars,
        figure_size_x=figure_size_x,
        figure_size_y=figure_size_y,
        x_plot_range_apparent=x_plot_range_apparent,
        y_plot_range_apparent=y_plot_range_apparent,
        x_plot_range_absolute=x_plot_range_absolute,
        y_plot_range_absolute=y_plot_range_absolute,
        isochrone_configuration_file=isochrone_configuration_file,
    )

    print(f'{Bcolors.OKGREEN}   Done{Bcolors.ENDC}')
