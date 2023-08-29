#! /usr/bin/python3
# -*- coding: utf-8 -*-

##################################################################################
#                               Script Parameters                                #
##################################################################################

#   Science spectrum file
science = 'star.FIT'

#   Directory of the dark frame for the stellar spectrum
darkframe_dir = 'darks/??s/'

#   Flatfield directory
flatfield_dir = 'flats/'

#   Directory of the dark frame for the flats
flatdark_dir = 'darks/??s/'

#   Region containing the science spectrum
specRegionStart = 495
specRegionEnd = 590

#   Sky background region (inside the slit)
bgSkyStart = 480
bgSkyEnd = 490

###
#   Plot range
#
#   Set the variables to '?' for an automatic resizing
lambdamin = '?'
lambdamax = '?'
# lambdamin = 3500.
# lambdamax = 5000.

###
#   Idents
#
#   Plot idents yes or no
plotident = 'yes'
# plotident = 'no'

#   File containing line identifications
#   lineFile   = ""
lineFile = "absorption_lines.dat"

##################################################################################
#          The following parameters usually do not need to be adjusted           #
##################################################################################

#   Image reduction mode
# mode = 'mean'
mode = 'median'

#   Calibration file
calib = 'calibration_spectrum.dat'

#   Data output
spectrumfile = 'stern_spectrum.pdf'
spectrumdata = 'stern_spectrum.dat'

masterdark_output_name = 'master_dark.fit'
masterflat_output_name = 'master_flat.fit'

#   Background region (outside the slit)
bgRegionStart = 0
bgRegionEnd = 0

##################################################################################
#                            Script Routines                                     #
##################################################################################

import os, pylab, time, sys, matplotlib
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from astropy.table import Table

matplotlib.rcParams['pdf.fonttype'] = 42


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


##################################################################################
#                               Functions                                        #
##################################################################################

print(bcolors.BOLD + "   Check input data" + bcolors.ENDC)

if not os.path.isfile(science):
    print(bcolors.FAIL + "   [ERROR] File containing spectrum doesn't exist. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if not os.path.exists(darkframe_dir):
    print(bcolors.FAIL + "   [ERROR] Direcory to Darkframes don't exist! " + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()
if os.path.isfile(darkframe_dir):
    print(bcolors.FAIL + "   [ERROR] Darkframes are linked to file, give the directory instead!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if not os.path.exists(flatfield_dir):
    print(bcolors.FAIL + "   [ERROR] Direcory to Flatfields don't exist! " + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()
if os.path.isfile(flatfield_dir):
    print(bcolors.FAIL + "   [ERROR] Flatfields are linked to file, give the directory instead!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if not os.path.exists(flatdark_dir):
    print(bcolors.FAIL + "   [ERROR] Direcory to Flatdarks don't exist! " + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()
if os.path.isfile(flatdark_dir):
    print(bcolors.FAIL + "   [ERROR] Flatdarks are linked to file, give the directory instead!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if not os.path.isfile(calib):
    print(bcolors.FAIL + "   [ERROR] Calibration file doesn't exist. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

flatdark_list = os.listdir(flatdark_dir)
flatfield_list = os.listdir(flatfield_dir)
darkframe_list = os.listdir(darkframe_dir)

if len(flatdark_list) == 0:
    print(bcolors.FAIL + "   [ERROR] Flatdark folder is empty. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if len(flatfield_list) == 0:
    print(bcolors.FAIL + "   [ERROR] Flatfield folder is empty. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if len(darkframe_list) == 0:
    print(bcolors.FAIL + "   [ERROR] Darkframe folder is empty. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

science_image = fits.open(science)
setup_dark_image = fits.open(darkframe_dir + "/" + darkframe_list[0])
setup_flat_image = fits.open(flatfield_dir + "/" + flatfield_list[0])
setup_flatdark_image = fits.open(flatdark_dir + "/" + flatdark_list[0])

science_image_data = science_image[0].data
size_science_x, size_science_y = len(science_image_data), len(science_image_data[1])

setup_dark_image_data = setup_dark_image[0].data
size_dark_x, size_dark_y = len(setup_dark_image_data), len(setup_dark_image_data[1])

setup_flat_image_data = setup_flat_image[0].data
size_flat_x, size_flat_y = len(setup_flat_image_data), len(setup_flat_image_data[1])

setup_flatdark_image_data = setup_flatdark_image[0].data
size_flatdark_x, size_flatdark_y = len(setup_flatdark_image_data), len(setup_flatdark_image_data[1])

setup_flatdark_image.close()
setup_dark_image.close()
setup_flat_image.close()

if size_flatdark_x != size_flat_x or size_flatdark_y != size_flat_y:
    print(bcolors.FAIL + "   [ERROR] Flatdarks and Flatfields don't have the same size, check!" + bcolors.ENDC)
    print(
        bcolors.FAIL + "   [ERROR] (" + size_flatdark_x + "x" + size_flatdark_y + ") and (" + size_flat_x + "x" + size_flat_y + ")" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if size_dark_x != size_flat_x or size_dark_y != size_flat_y:
    print(bcolors.FAIL + "   [ERROR] Flatdarks and Flatfields don't have the same size, check!" + bcolors.ENDC)
    print(
        bcolors.FAIL + "   [ERROR] (" + size_flatdark_x + "x" + size_flatdark_y + ") and (" + size_flat_x + "x" + size_flat_y + ")" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if size_dark_x != size_science_x or size_dark_y != size_science_y:
    print(bcolors.FAIL + "   [ERROR] Flatdarks and Flatfields don't have the same size, check!" + bcolors.ENDC)
    print(
        bcolors.FAIL + "   [ERROR] (" + size_flatdark_x + "x" + size_flatdark_y + ") and (" + size_flat_x + "x" + size_flat_y + ")" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

### Data reduction ###

print(bcolors.BOLD + "   Create Masterdark from " + str(len(darkframe_list)) + " files" + bcolors.ENDC)
dark_array = np.zeros((size_dark_x, size_dark_y, len(darkframe_list)))
for element in darkframe_list:
    dark_image = fits.open(darkframe_dir + "/" + element)
    dark_image_data = dark_image[0].data
    darkHeader = dark_image[0].header
    dark_array[..., darkframe_list.index(element)] = dark_image_data
    dark_image.close()
if mode == "median":
    dark_out_array = np.median(dark_array.astype(int), axis=2)
if mode == "mean":
    dark_out_array = np.mean(dark_array.astype(int), axis=2)

hdu = fits.PrimaryHDU(dark_out_array, darkHeader)
hdu.scale('int16', '', bzero=0)
hduList = fits.HDUList([hdu])
hduList.writeto(masterdark_output_name, output_verify='exception', overwrite=True)

print(bcolors.BOLD + "   Create Flatdark from " + str(len(flatdark_list)) + " files" + bcolors.ENDC)
flatdark_array = np.zeros((size_flat_x, size_flat_y, len(flatdark_list)))
for element in flatdark_list:
    flatdark_image = fits.open(flatdark_dir + "/" + element)
    flatdark_image_data = flatdark_image[0].data
    flatdark_array[..., flatdark_list.index(element)] = flatdark_image_data
    flatdark_image.close()
if mode == "mean":
    flatdark_out_array = np.mean(flatdark_array.astype(int), axis=2)
if mode == "median":
    flatdark_out_array = np.median(flatdark_array.astype(int), axis=2)

print(bcolors.BOLD + "   Create Masterflat from " + str(len(flatfield_list)) + " files" + bcolors.ENDC)
flatfield_array = np.zeros((size_flat_x, size_flat_y, len(flatfield_list)))
for element in flatfield_list:
    flatfield_image = fits.open(flatfield_dir + "/" + element)
    flatfield_image_data = flatfield_image[0].data
    flatHeader = flatfield_image[0].header
    flatfield_array[..., flatfield_list.index(element)] = flatfield_image_data - flatdark_out_array
    flatfield_image.close()
if mode == "mean":
    flatfield_out_array = np.mean(flatfield_array.astype(int), axis=2)
if mode == "median":
    flatfield_out_array = np.median(flatfield_array.astype(int), axis=2)

hdu = fits.PrimaryHDU(flatfield_out_array, flatHeader)
hdu.scale('int16', '', bzero=0)
hduList = fits.HDUList([hdu])
hduList.writeto(masterflat_output_name, overwrite=True)

print(bcolors.BOLD + "   Apply Drakframe correction to science spectrum" + bcolors.ENDC)
science_array = np.asarray(science_image_data)
sience_reduced_array = np.zeros((size_flat_x, size_flat_y, len(flatdark_list)))
sience_reduced_array = (science_array - dark_out_array)  # / (flatfield_out_array/np.amax(flatfield_out_array ))
# sience_reduced_array = science_array
# flatfield_out_array = flatfield_out_array/np.amax(flatfield_out_array)

print(bcolors.BOLD + "   Extract spectral data" + bcolors.ENDC)
# spectrum star
spec = []
# spectrum flatfield
spec2 = []
# Sky background spectrum
skyspec = []

for i in range(0, len(sience_reduced_array[0])):
    spec_tmp, bg_tmp = 0, 0
    for j in range(specRegionStart, specRegionEnd):
        spec_tmp += sience_reduced_array[j][i]
    # for k in range(bgRegionStart,bgRegionEnd):
    # bg_tmp += sience_reduced_array[k][i]
    spec_tmp /= abs((specRegionEnd - specRegionStart))
    # bg_tmp /=  abs((bgRegionEnd-bgRegionStart))
    spec.append(spec_tmp - bg_tmp)

for i in range(0, len(flatfield_out_array[0])):
    spec_tmp, bg_tmp = 0, 0
    for j in range(specRegionStart, specRegionEnd):
        spec_tmp += flatfield_out_array[j][i]
    # for k in range(bgRegionStart,bgRegionEnd):
    # bg_tmp += flatfield_out_array[k][i]
    spec_tmp /= abs((specRegionEnd - specRegionStart))
    # bg_tmp /=  abs((bgRegionEnd-bgRegionStart))
    spec2.append(spec_tmp - bg_tmp)

for i in range(0, len(sience_reduced_array[0])):
    skybg, bg_tmp = 0, 0
    for j in range(bgSkyStart, bgSkyEnd):
        skybg += sience_reduced_array[j][i]
    # for k in range(bgRegionStart,bgRegionEnd):
    # bg_tmp += sience_reduced_array[k][i]
    skybg /= abs((bgSkyStart - bgSkyEnd))
    # bg_tmp /=  abs((bgRegionEnd-bgRegionStart))
    skyspec.append(skybg - bg_tmp)

print(bcolors.BOLD + "   Apply Flatfield correction to science spectrum" + bcolors.ENDC)

# normalize the flatfield
spec2 = np.asarray(spec2) / max(spec2)
# apply flatfield to science spectrum and sky spectrum
spec = spec / spec2
skyspec = skyspec / spec2

print(bcolors.BOLD + "   Perform wavelength calibration" + bcolors.ENDC)

calibfile = open(calib, 'r')
wavelengthrange = []
for line in calibfile:
    liste = line.split()
    if len(liste) == 0:
        continue
    if lambdamin != '?' and float(liste[0]) < lambdamin:
        spec = spec[1:]
        spec2 = spec2[1:]
        continue
    if lambdamax != '?' and float(liste[0]) > lambdamax:
        spec = spec[:-1]
        spec2 = spec2[:-1]
        continue
    wavelengthrange.append(float(liste[0]))

print(bcolors.BOLD + "   Create spectral plot " + bcolors.OKBLUE + spectrumfile + bcolors.ENDC)

### Plotting the spectrum ###
fig1 = plt.figure(figsize=(20, 10))

font = {'family': 'serif',
        'color': 'red',
        'weight': 'normal',
        'size': 15,
        }

# Setting plot labels
plt.xlabel(r'$\lambda\,[\AA]$')
plt.ylabel('Relative flux')

# Setting plot ranges
yoffset = (max(spec) - min(spec)) * 0.05
pylab.ylim([min(spec) - yoffset, max(spec) + yoffset])

plotoffset = (float(max(wavelengthrange)) - float(min(wavelengthrange))) * 0.01
pylab.xlim(min(wavelengthrange) - plotoffset, max(wavelengthrange) + plotoffset)

# Plot the actual data
plt.plot(wavelengthrange, spec - skyspec, 'b-')
# plt.plot(wavelengthrange,spec2,'r-')
# plt.show()

# Setting plotpositions for ident lines
if plotident == 'yes' and lineFile != "":
    plotminimum = min(spec) - yoffset
    plotmaximum = max(spec) + yoffset
    plotheigth = plotmaximum - plotminimum
    plotmiddleplot = (plotminimum + plotmaximum) / 2.0
    plotupperplot = plotminimum + 0.80 * plotheigth
    plotlowerplot = plotminimum + 0.20 * plotheigth
    plotuppercut1 = plotminimum + 0.70 * plotheigth
    plotlowercut1 = plotminimum + 0.30 * plotheigth
    plotuppercut2 = plotminimum + 0.68 * plotheigth
    plotlowercut2 = plotminimum + 0.32 * plotheigth

    # interpolate on data to find point for ident
    f2 = interp1d(wavelengthrange, spec)
    lines = open(lineFile, "r")

    # plot idents to figure
    for line in lines:
        liste = line.split()
        if len(liste) == 1:
            print(
                bcolors.WARNING + "     [WARNING] Broken identification found as '" + line + "', must consist of [wavelenght(s) + name]. I will skip this one." + bcolors.ENDC)
            continue
        try:
            float(liste[0])
        except ValueError:
            print(
                bcolors.WARNING + "     [WARNING] Broken identification found as '" + line + "', first entry not a number. I will skip this one." + bcolors.ENDC)
            continue

        #   Single ident plot
        if len(liste) == 2:
            ident_line = liste
            wave_line = float(ident_line[0])
            if float(wave_line) >= min(wavelengthrange) and float(wave_line) <= max(wavelengthrange):
                if f2(wave_line) <= plotmiddleplot:
                    plt.plot(
                        [wave_line, wave_line],
                        [plotupperplot, f2(wave_line)],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        wave_line,
                        plotupperplot,
                        ident_line[1],
                        rotation=90,
                        ha='center',
                        va='bottom',
                        fontdict=font,
                    )
                else:
                    plt.plot(
                        [wave_line, wave_line],
                        [f2(wave_line), plotlowerplot],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        wave_line,
                        plotlowerplot,
                        ident_line[1],
                        rotation=90,
                        ha='center',
                        va='top',
                        fontdict=font,
                    )
        #   Multi ident plot
        if len(liste) > 2:
            points = []
            for i in liste[:-1]:
                points.append(float(i))
            pointcenter = sum(points) / float(len(points))
            ident_name = str(liste[-1:])
            # print(points," give ",pointcenter," bei ",liste[-1:])
            if max(points) <= max(wavelengthrange) and min(points) >= min(wavelengthrange):
                if f2(pointcenter) <= plotmiddleplot:
                    plt.plot(
                        [pointcenter, pointcenter],
                        [plotupperplot, plotuppercut1],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        pointcenter,
                        plotupperplot,
                        ident_name[2:-2],
                        rotation=90,
                        ha='center',
                        va='bottom',
                        fontdict=font,
                    )
                    for element in points:
                        plt.plot(
                            [element, element],
                            [f2(element), plotuppercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.5,
                        )
                        plt.plot(
                            [pointcenter, element],
                            [plotuppercut1, plotuppercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.5,
                        )
                if f2(pointcenter) > plotmiddleplot:
                    plt.plot(
                        [pointcenter, pointcenter],
                        [plotlowerplot, plotlowercut1],
                        color='r',
                        linestyle='-',
                        linewidth=1.5,
                    )
                    plt.text(
                        pointcenter,
                        plotlowerplot,
                        ident_name[2:-2],
                        rotation=90,
                        ha='center',
                        va='top',
                        fontdict=font,
                    )
                    for element in points:
                        plt.plot(
                            [element, element],
                            [f2(element), plotlowercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.5,
                        )
                        plt.plot(
                            [pointcenter, element],
                            [plotlowercut1, plotlowercut2],
                            color='r',
                            linestyle='-',
                            linewidth=1.5,
                        )
    lines.close()

# Write the plot file

plt.savefig(spectrumfile, bbox_inches='tight')

plt.clf()

### Plotting the flatfield ###
fig2 = plt.figure(figsize=(20, 10))
# Setting plot labels
plt.xlabel(r'$\lambda\,[\AA]$')
plt.ylabel('Relative flux')

# Setting plot ranges
yoffset = (max(spec2) - min(spec2)) * 0.05
pylab.ylim([min(spec2) - yoffset, max(spec2) + yoffset])

plotoffset = (float(max(wavelengthrange)) - float(min(wavelengthrange))) * 0.01
pylab.xlim(min(wavelengthrange) - plotoffset, max(wavelengthrange) + plotoffset)

# Plot the actual data
plt.plot(wavelengthrange, spec2, 'r-')

print(bcolors.BOLD + "   Create flatfield plot " + bcolors.OKBLUE + 'flatfield.pdf' + bcolors.ENDC)

# Write the plot file
plt.savefig('flatfield.pdf', bbox_inches='tight')

print(bcolors.BOLD + "   Write spectrum to file " + bcolors.OKBLUE + spectrumdata + bcolors.ENDC)

#   Write spectrum table
tbl = Table(names=['wave', 'flux', ], data=[wavelengthrange, spec, ])
tbl.write(spectrumdata, format='ascii', overwrite=True)

# os.system('touch '+spectrumdata)
# specdatafile = open(spectrumdata,'w')

# for i in range(len(wavelengthrange)):
# specdatafile.write( str(wavelengthrange[i]) + "   " +str(spec[i]))

print(bcolors.OKGREEN + "   Done" + bcolors.ENDC)
