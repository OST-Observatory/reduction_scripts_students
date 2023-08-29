#! /usr/bin/python3
# -*- coding: utf-8 -*-

##################################################################################
#                               Script Parameters                                #
##################################################################################

#   Name of files
spectrum_1_file_name = ""
spectrum_2_file_name = ""

#   Start row and end row in the first image with the name "spectrum_1_file_name"
row_1_start = 0
row_1_end = 1

#   Start row and end row in the second image with the name "spectrum_2_file_name"
row_2_start = 0
row_2_end = 1

#   Offset of the spectra along the x-axis (in pixel)
offset_spectrum_1 = 0
offset_spectrum_2 = 0

#   Shift of the spectra along the y-axis (in flux units) shift is applied after normalization
shift_spectrum_1 = 0.0
shift_spectrum_2 = 0.0

#   Normalize the spectra? well, kind of...
# normalize = "no"
normalize = "yes"

#   Zoom into a region along x axis?
# zoom = 'yes'
zoom = 'no'
# zoom range
x_start = 2000
x_end = 2200

#   Name of output file
filename = "spectrum"  # default spectrum
filetype = "pdf"  # supported filetypes: png, pdf, ps, eps and svg, default: pdf

#   Write output spectrum in two separate files?
two_files = 'yes'
# two_files = 'no'
# in case of one file, give file name:
output_file_name = 'spectrum.dat'

##################################################################################
#                                Script Routines                                 #
##################################################################################

import os, sys, string, pylab
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


#   Define font colors for terminal output
class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


#   Import data from files
print(Bcolors.BOLD + "   Check input parameter" + Bcolors.ENDC)

#   Check if spectrum files are accessible
if not os.path.isfile(spectrum_1_file_name):
    print(Bcolors.FAIL + "   [ERROR] I can't find " + spectrum_1_file_name + ", check file name!" + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()

if not os.path.isfile(spectrum_2_file_name):
    print(Bcolors.FAIL + "   [ERROR] I can't find " + spectrum_2_file_name + ", check file name!" + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()

#   Check if offsets and are numbers
try:
    float(offset_spectrum_1)
except ValueError:
    print(Bcolors.FAIL + "   [ERROR] x offset for first spectrum is not a number" + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()
try:
    float(offset_spectrum_2)
except ValueError:
    print(Bcolors.FAIL + "   [ERROR] x offset for second spectrum is not a number" + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()
try:
    float(shift_spectrum_1)
except ValueError:
    print(Bcolors.FAIL + "   [ERROR] y shift for first spectrum is not a number" + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()
try:
    float(shift_spectrum_2)
except ValueError:
    print(Bcolors.FAIL + "   [ERROR] y shift for second spectrum is not a number" + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()

if (type(row_1_start).__name__ != "int" or
        type(row_2_start).__name__ != "int" or
        type(row_1_end).__name__ != "int" or
        type(row_2_end).__name__ != "int"):
    print(Bcolors.FAIL + "   [ERROR] rows need to be integers " + Bcolors.ENDC)
    print(Bcolors.FAIL + "   [ERROR] ABORT" + Bcolors.ENDC)
    sys.exit()

print(Bcolors.BOLD + "   Read Fits files" + Bcolors.ENDC)
#   Open fits files and read the data sections into 2-dimensional arrays
#   File with spectrum
spectrum_1_file = fits.open('%s' % spectrum_1_file_name)
spectrum_1_data = spectrum_1_file[0].data

#   Dark frame for the spectrum
spectrum_2_file = fits.open('%s' % spectrum_2_file_name)
spectrum_2_data = spectrum_2_file[0].data

#   Get dimensions of the data
n_data_columns = len(spectrum_1_data[0])  # number of columns in the file
n_data_lines = len(spectrum_1_data)  # number of lines in the file

#   Read lines that contain spectral information
print(Bcolors.BOLD + "   Extract spectra " + Bcolors.ENDC)
extracted_spectrum_1 = list(
    np.average(spectrum_1_data[row_1_start:row_1_end], axis=0)
)
extracted_spectrum_2 = list(
    np.average(spectrum_2_data[row_2_start:row_2_end], axis=0)
)
# spec1 = []
# tmp = 0
# for j in range(0, dataColumns):
# for i in range(row1start,row1end):
# tmp += spec1FileData[i][j]
# tmp /= (row1end-row1start)
# spec1.append(tmp)

# spec2 = []
# tmp = 0
# for j in range(0, dataColumns):
# for i in range(row2start,row2end):
# tmp = tmp + spec2FileData[i][j]
# tmp = tmp / (row2end-row2start)
# spec2.append(tmp)

# print(spec1,spec2)


#   Set automatic plot range
if len(extracted_spectrum_1) + offset_spectrum_1 > len(extracted_spectrum_2) + offset_spectrum_2:
    plot_x_max = 1.02 * (len(extracted_spectrum_1) + offset_spectrum_1)
else:
    plot_x_max = 1.02 * (len(extracted_spectrum_2) + offset_spectrum_2)

if offset_spectrum_1 < offset_spectrum_2:
    plot_x_min = offset_spectrum_1 - .02 * (len(extracted_spectrum_1) + offset_spectrum_2)
else:
    plot_x_min = offset_spectrum_2 - .02 * (len(extracted_spectrum_1) + offset_spectrum_1)

plt.figure(figsize=(10, 5))

print(Bcolors.BOLD + "   Plot spectra " + Bcolors.ENDC)

#   Adjust plot range if zoom is selected
if zoom == 'yes':
    try:
        float(x_start)
        float(x_end)
    except ValueError:
        print(Bcolors.WARNING + "   [WARNING] Zoom requested but range is not a number" + Bcolors.ENDC)
        print(Bcolors.WARNING + "             -> automatic plot range used" + Bcolors.ENDC)
    else:
        if x_start > x_end:
            a = x_end
            x_end = x_start
            x_start = a
        plot_x_min, plot_x_max = float(x_start), float(x_end)
        print(
            f"{Bcolors.BOLD}   Zoom in on range ({Bcolors.OKBLUE}{x_start},{x_end}){Bcolors.ENDC}"
        )

plt.xlim([plot_x_min, plot_x_max])

#   Apply the normalization by dividing by max value amd plot spectrum
if normalize == "yes":
    if shift_spectrum_1 > shift_spectrum_2:
        pylab.ylim([0, 1.1 + shift_spectrum_1])
    else:
        pylab.ylim([0, 1.1 + shift_spectrum_2])
    plt.plot(
        np.arange(offset_spectrum_1, len(extracted_spectrum_1) + offset_spectrum_1),
        extracted_spectrum_1 / max(extracted_spectrum_1) + shift_spectrum_1,
        'r-',
        linewidth=.5,
        label=os.path.split(spectrum_1_file_name)[-1][:-4],
    )
    plt.plot(
        np.arange(offset_spectrum_2, len(extracted_spectrum_2) + offset_spectrum_2),
        extracted_spectrum_2 / max(extracted_spectrum_2) + shift_spectrum_2,
        'b-',
        linewidth=.5,
        label=os.path.split(spectrum_2_file_name)[-1][:-4],
    )
    extracted_spectrum_1 = extracted_spectrum_1 / max(extracted_spectrum_1)
    extracted_spectrum_2 = extracted_spectrum_2 / max(extracted_spectrum_2)
else:
    plt.plot(
        np.arange(offset_spectrum_1, len(extracted_spectrum_1) + offset_spectrum_1),
        extracted_spectrum_1 + shift_spectrum_1,
        'r-',
        linewidth=.5,
    )
    plt.plot(
        np.arange(offset_spectrum_2, len(extracted_spectrum_2) + offset_spectrum_2),
        extracted_spectrum_2 + shift_spectrum_2,
        'b-',
        linewidth=.5,
    )

#   Set output filenames
if filename == "" or filename == "?":
    print(Bcolors.WARNING + "   [WARNING] No filename given, use default 'spectrum'" + Bcolors.ENDC)
    filename = "spectrum"

if filetype != "pdf" and filetype != "png" and filetype != "eps" and filetype != "ps" and filetype != "svg":
    print(Bcolors.WARNING + "   [WARNING] No or unknown filetype given, use default 'pdf'" + Bcolors.ENDC)
    filetype = "pdf"

#   Plot legend
plt.legend(loc=3, ncol=1)
# adjust ticks to appear all around
plt.tick_params(top=True, right=True, which='both', direction='in')
plt.minorticks_on()

# os.system('rm tmpdata')
# set plot labels
plt.xlabel('Pixel')
if normalize == 'yes':
    plt.ylabel('Normalized intensity')
else:
    plt.ylabel('Intensity')

print(Bcolors.BOLD + "   Write " + Bcolors.OKBLUE + filename + '.' + filetype + Bcolors.ENDC)
#   Save figure
plt.savefig(filename + '.' + filetype, format=filetype, bbox_inches='tight')

#   Write spectral data into data set
if two_files == 'no':
    print(Bcolors.BOLD + '   Write output/' + Bcolors.OKBLUE + os.path.split(spectrum_1_file_name)[-1][
                                                               :-3] + "dat" + Bcolors.ENDC)
    os.system("mkdir -p output")
    os.system('touch output/' + os.path.split(spectrum_1_file_name)[-1][:-3] + "dat")
    f = open('output/' + os.path.split(spectrum_1_file_name)[-1][:-3] + "dat", 'w')
    for i in range(0, len(extracted_spectrum_1)):
        j = i + offset_spectrum_1
        f.write('%i\t%f\n' % (j, extracted_spectrum_1[i]))
    f.close()

    print(Bcolors.BOLD + '   Write output/' + Bcolors.OKBLUE + os.path.split(spectrum_2_file_name)[-1][
                                                               :-3] + "dat" + Bcolors.ENDC)
    os.system('touch output/' + os.path.split(spectrum_2_file_name)[-1][:-3] + "dat")
    f = open('output/' + os.path.split(spectrum_2_file_name)[-1][:-3] + "dat", 'w')
    for i in range(0, len(extracted_spectrum_2)):
        j = i + offset_spectrum_2
        f.write('%i\t%f\n' % (j, extracted_spectrum_2[i]))
    f.close()

else:
    print(Bcolors.BOLD + '   Write output/' + Bcolors.OKBLUE + output_file_name + Bcolors.ENDC)
    os.system("mkdir -p output")
    os.system('touch output/' + output_file_name)
    f = open('output/' + output_file_name, 'w')
    for i in range(0, len(extracted_spectrum_2)):
        j = i + offset_spectrum_2
        f.write(f'{j:>5} {extracted_spectrum_1[i]:>25} {str(extracted_spectrum_2[i]):>25}\n')
    f.close()

print(Bcolors.OKGREEN + "   Done" + Bcolors.ENDC)
