#! /usr/bin/python3
# -*- coding: utf-8 -*-

############################################################################
#             Configuration: modify the file in this section               #
############################################################################

#   Name of the file with the wavelength calibration spectrum
calibFileName = "master_wave.fit"

#   Region (rows on the image) containing the calibration spectrum
specRegionStart = 476
specRegionEnd = 487

#   Background region (rows on the image) -> need to be outside of the slits
bgRegionStart = 815
bgRegionEnd = 870

############################################################################
#                Additional options: only edit if necessary                #
############################################################################

#   Numerical parameters for line fit
#   (width of emission lines in pixel)
emissionlinewidth = 1.5
#   Number of emission lines to look for
emission_lines_to_find = 80

#   Names of output files
waveFile = "calibration_spectrum.dat"
outputwave = "calibration_fit.pdf"
outputselect = "calibration_selection.pdf"

#   Expected lines in calibration spectrum
#   Hg & Ar lamps
linelist = [
    3650.158,
    4046.56,
    4077.837,
    4358.33,
    5460.75,
    5769.610,
    5790.670,
    6965.431,
    7067.218,
    7147.042,
    7272.936,
    7383.980,
    7503.86,
    7635.44,
    7724.20,
    7948.176,
    8006.157,
    8103.693,
    8115.81,
    8264.522,
    8408.210,
    8521.442,
    8667.944,
    # 9122.967,
]

#   Ar & Ne lamps
linelist = [
#     4131.72
#     4158.59,
    4200.67,
    4259.36,
    4277.53,
    # 4348.06,
    4510.73,
    # 4545.05,
    4609.57,
    4657.90,
    4764.87,
    4861.35,
    4879.86,
    5852.49,
    6143.06,
    6402.25,
    6562.79,
    6677.28,
    6965.43,
    7067.22,
    7272.94,
    7383.98,
    7503.87,
    7635.11,
    7723.76,
    7948.18,
    8014.79,
    8103.69,
    8115.31,
    8264.52,
    8408.21,
    8418.43,
    8521.44,
    9122.97,
]


############################################################################
#                               Libraries                                  #
############################################################################

import os, pylab, time, sys, matplotlib
from numpy import matrix
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from decimal import Decimal

############################################################################
#                           Routines & definitions                         #
############################################################################

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


def gaus(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def click_point(event):
    if event.button == 1:
        # print('left click on '+str(event.xdata),str(event.ydata))
        # fig1.suptitle('last click: left')

        if len(linelist_tmp) > 0:
            # print('ECHO')
            # print(linelist_tmp)
            distance = []
            for i in range(len(maxindexlist)):
                distance.append(np.sqrt(((event.xdata - maxindexlist[i]) / len(spec)) ** 2 + (
                            (event.ydata - maxvaluelist[i]) / ((max(spec) + 0.1 * max(spec)) * (16 / 9))) ** 2))
                # distance.append(np.sqrt((event.xdata)**2+(event.ydata)*ratio)**2)
            indexx = distance.index(min(distance))
            # print((((event.xdata-maxindexlist[indexx])),((event.ydata-maxvaluelist[indexx]))))
            # print((event.xdata-maxindexlist[indexx])/len(spec),(event.ydata-maxvaluelist[indexx])/(max(spec)+0.1*max(spec)),max(spec),len(spec))
            # print(min(distance),indexx)
            fit_y.append(linelist_tmp[0])
            fit_x.append(maxindexlist[indexx])

            plt.text(
                maxindexlist[indexx] + len(spec) / 100,
                maxvaluelist[indexx],
                str(linelist_tmp[0]),
                va="center",
                ha="left",
            )
            linelist_tmp.pop(0)

            if len(linelist_tmp) != 0:
                fig1.suptitle(
                    r'$\lambda$ = ' + str(linelist_tmp[0]) + r' $\AA$',
                    color='r',
                    fontsize=20,
                )
            else:
                fig1.suptitle(
                    'All lines done',
                    color='r',
                    fontsize=20,
                )
            plt.plot(
                maxindexlist[indexx],
                maxvaluelist[indexx],
                'bo',
            )

            # plt.plot([maxindexlist[indexx],event.xdata],[maxvaluelist[indexx],event.ydata],'b-')
            event.canvas.draw()

            maxindexlist.pop(indexx)
            maxvaluelist.pop(indexx)

            if len(linelist_tmp) == 0:
                time.sleep(1.0)
                fig1.suptitle("")
                plt.savefig(outputselect, bbox_inches="tight")
                plt.close()

        else:
            fig1.suptitle(
                'All lines done',
                color='r',
                fontsize=20,
            )
            time.sleep(1.0)
            fig1.suptitle("")
            plt.savefig(outputselect, bbox_inches="tight")
            plt.close()

    if event.button == 3:
        # print('right click on '+str(event.xdata),str(event.ydata))
        if len(linelist_tmp) > 0:
            linelist_tmp.pop(0)
            # fig1.title('asd',color='g')

            if len(linelist_tmp) != 0:
                fig1.suptitle(
                    r'$\lambda$ = ' + str(linelist_tmp[0]) + r' $\AA$',
                    color='r',
                    fontsize=20,
                )
            else:
                fig1.suptitle(
                    'All lines done',
                    color='r',
                    fontsize=20,
                )
            event.canvas.draw()

            if len(linelist_tmp) == 0:
                time.sleep(1.0)
                fig1.suptitle("")
                plt.savefig(outputselect, bbox_inches="tight")
                plt.close()
        else:
            fig1.suptitle(
                'All lines done',
                color='r',
                fontsize=20,
            )
            time.sleep(1.0)
            fig1.suptitle("")
            plt.savefig(outputselect, bbox_inches="tight")
            plt.close()


def press_button(event):
    if event.key == 'q':
        fig1.suptitle("")
        plt.savefig(outputselect, bbox_inches="tight")
        plt.close()


def poly3(x):
    return m1 * (x ** 3) + m2 * (x ** 2) + m3 * x + b


############################################################################
#                                  Main                                    #
############################################################################

print(bcolors.BOLD + "   Check input data" + bcolors.ENDC)

if os.path.isfile(calibFileName) == False:
    print(bcolors.FAIL + "   [ERROR] Calibration files doesn't exist. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

try:
    int(emission_lines_to_find)
except ValueError:
    print(bcolors.FAIL + "   [ERROR] Number of emission lines to find is not a number. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

print(bcolors.BOLD + "   Read fits file" + bcolors.ENDC)

calibFile = fits.open(calibFileName)
calibFileData = calibFile[0].data

## get dimensions of the data
calibFileDataColumns = len(calibFileData[0])  # number of columns in the file
calibFileDataLines = len(calibFileData)  # number of lines in the file
N = calibFileDataColumns

if bgRegionStart > calibFileDataLines or bgRegionEnd > calibFileDataLines:
    print(bcolors.FAIL + "   [ERROR] Background region outside of image. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

if specRegionStart > calibFileDataLines or specRegionEnd > calibFileDataLines:
    print(bcolors.FAIL + "   [ERROR] Background region outside of image. Check!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

print(bcolors.BOLD + "   Extract spectral data" + bcolors.ENDC)
linelist_tmp = []
linelist_tmp = linelist
spec = []
for i in range(0, N):
    spec_tmp, bg_tmp = 0, 0
    for j in range(specRegionStart, specRegionEnd):
        spec_tmp += calibFileData[j][i]
    for k in range(bgRegionStart, bgRegionEnd):
        bg_tmp += calibFileData[k][i]
    spec_tmp /= (specRegionEnd - specRegionStart)
    bg_tmp /= (bgRegionEnd - bgRegionStart)
    spec.append(spec_tmp - bg_tmp)

specrange = np.arange(len(spec))

counter = 0
spec_tmp = spec.copy()
asd = np.asarray(spec)
fig1 = plt.figure(num='Select line for the given wavelength', figsize=(16, 9))  # figsize=(80, 60))

plt.xlabel('Pixel')
plt.ylabel('Counts')
plt.plot(specrange, spec, "b-", linewidth=0.3)
plt.xlim([0, len(spec)])
plt.ylim([0, max(spec) + 0.1 * max(spec)])

maxindexlist = []
maxvaluelist = []
distance = []

print(bcolors.BOLD + "   Search for " + str(emission_lines_to_find) + " emissions lines" + bcolors.ENDC)

for i in range(0, emission_lines_to_find):
    max_value = max(spec_tmp)
    max_index = spec_tmp.index(max_value)
    plt.plot(max_index, max_value, "ro", linewidth=5)
    for j in range(0, len(spec_tmp)):
        spec_tmp[j] = spec_tmp[j] - gaus(j, max_value, max_index, emissionlinewidth)
    maxindexlist.append(max_index)
    maxvaluelist.append(max_value)

ratio = (max(maxindexlist) - min(maxindexlist)) / (max(maxvaluelist) - min(maxvaluelist))
fit_x = []
fit_y = []

fig1.suptitle(r'$\lambda$ = ' + str(linelist_tmp[0]) + r' $\AA$', color='r', fontsize=20)

fig1.canvas.mpl_connect('button_press_event', click_point)
fig1.canvas.mpl_connect('key_press_event', press_button)

# print("I will open a plot window of the extracted spectrum now. The lines i was able to find are marked, maybe not all of them are real. Left-Click on the marker for the wavelength I will display. Right-click to skip a wavelength.")

print(bcolors.BOLD + "   Emissions lines identified, calibration plot opened" + bcolors.ENDC)

print(bcolors.OKBLUE + "   [ACTION REQUIRED] Select matching emission line with displayed wavelength" + bcolors.ENDC)
print(bcolors.OKBLUE + "   [ACTION REQUIRED] Left Click: select line" + bcolors.ENDC)
print(bcolors.OKBLUE + "   [ACTION REQUIRED] Right Click: skip current wavelength" + bcolors.ENDC)
print(bcolors.OKBLUE + "   [ACTION REQUIRED] Select at least 4 lines" + bcolors.ENDC)
print(bcolors.OKBLUE + "   [ACTION REQUIRED] Press 'q' or close plot window to finish selection" + bcolors.ENDC)

plt.show()
plt.clf()
plt.close('all')

fig2 = plt.figure(num='Wavelength calibration fit')

plt.ylabel(r'$\lambda\,[\AA]$')
plt.xlabel('Pixel')
# def press_key(event):
# if event.key == 'q':
# plt.close()

if len(fit_x) < 4:
    print(bcolors.FAIL + "   [ERROR] Not enough lines selected for a fit. Try again!" + bcolors.ENDC)
    print(bcolors.FAIL + "   [ERROR] ABORT" + bcolors.ENDC)
    sys.exit()

print(bcolors.BOLD + "   Fit selected lines (" + str(len(fit_x)) + ") with 3rd order polynomial " + bcolors.ENDC)
m1, m2, m3, b = pylab.polyfit(fit_x, fit_y, 3)
# print(m1,m2,m3,b)
print(bcolors.OKGREEN + "   Fit: " + 'y = ' + str('%.2E' % Decimal(m1)) + " x^3$ + " + str(
    '%.2E' % Decimal(m2)) + " x^2$ + " + str('%.2E' % Decimal(m3)) + "x + " + str('{:7.2f}'.format(b)) + bcolors.ENDC)

fig2.suptitle('$\lambda$ = ' + str('%.2E' % Decimal(m1)) + r" $\times$ pxl$^3$ + " + str(
    '%.2E' % Decimal(m2)) + r" $\times$ pxl$^2$ + " + str('%.2E' % Decimal(m3)) + r" $\times$ pxl + " + str(
    '{:7.2f}'.format(b)), fontsize=15)
# fig2.canvas.mpl_connect('key_press_event',press_key)

plt.plot(specrange, poly3(specrange), 'r--', linewidth=0.5)
plt.plot(fit_x, fit_y, 'go')

print(bcolors.BOLD + "   Fit is displayed in external window" + bcolors.ENDC)
print(bcolors.BOLD + "   This plot will be saved as " + bcolors.OKBLUE + outputwave + bcolors.ENDC)
print(bcolors.OKBLUE + "   [ACTION REQUIRED] Close window to continue" + bcolors.ENDC)

plt.show(block='Wavelength calibration fit')

print(bcolors.BOLD + "   Write calibrated calibration spectrum to " + bcolors.OKBLUE + waveFile + bcolors.ENDC)
plt.clf()
plt.close()

plt.ylabel(r'$\lambda\,[\AA]$')
plt.xlabel('Pixel')
plt.plot(specrange, poly3(specrange), 'r--', linewidth=0.5)
plt.plot(fit_x, fit_y, 'go')
plt.savefig(outputwave)

os.system('touch ' + waveFile)
specout = open(waveFile, 'w')

for i in range(len(specrange)):
    specout.write(str(poly3(specrange[i])) + '  ' + str(asd[i]) + "\n")
print(bcolors.OKGREEN + "   DONE" + bcolors.ENDC)
# print(fit_x,fit_y)
# print(b,m)
