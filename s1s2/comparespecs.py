#! /usr/bin/python3
# -*- coding: utf-8 -*-

##################################################################################
############################ Script Parameters ###################################
##################################################################################

# name of files
spec1FileName     = ""
spec2FileName     = ""


# start row and end row in the first image with the name "spec1FileName"
row1start = 0
row1end   = 1

# start row and end row in the second image with the name "spec2FileName"
row2start = 0
row2end   = 1

# offst of the spectra along the x-axis (in pixel)
offset_spec1 = 0
offset_spec2 = 0

# shift of the spectra along the y-axis (in flux units) shift is applied after normalization
shift_spec1 = 0.0
shift_spec2 = 0.0

# normalize the spectra? well, kind of...
#normalize = "no"
normalize = "yes"

#zoom into a region along x axis?
#zoom = 'yes'
zoom = 'no'
#zoom range
x_start = 2000
x_end   = 2200

# name of output file
filename  = "spectrum" #default spectrum
filetype  = "pdf" #supported filetypes: png, pdf, ps, eps and svg, default: pdf
#write outputspectrum in two separate files?
two_files = 'yes'
#two_files = 'no'
#in case of one file, give file name:
output_file_name = 'spectrum.dat'

##################################################################################
############################ Script Routines ####################################
##################################################################################

import os, sys, string,pylab
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

#Define font colors for terminal output
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

# import data from files
print(bcolors.BOLD+"   Check input parameter"+bcolors.ENDC)

#check if spectrum files are accessable
if not os.path.isfile(spec1FileName):
    print(bcolors.FAIL+"   [ERROR] I can't find "+spec1FileName+", check file name!"+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()
    
if not os.path.isfile(spec2FileName):
    print(bcolors.FAIL+"   [ERROR] I can't find "+spec2FileName+", check file name!"+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()

#check if offsets and are numbers
try: 
    float(offset_spec1)
except ValueError:
    print(bcolors.FAIL+"   [ERROR] x offset for first spectrum is not a number"+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()
try: 
    float(offset_spec2)
except ValueError:
    print(bcolors.FAIL+"   [ERROR] x offset for second spectrum is not a number"+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()
try: 
    float(shift_spec1)
except ValueError:
    print(bcolors.FAIL+"   [ERROR] y shift for first spectrum is not a number"+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()
try: 
    float(shift_spec2)
except ValueError:
    print(bcolors.FAIL+"   [ERROR] y shift for second spectrum is not a number"+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit()  
    
if type(row1start).__name__ != "int" or type(row2start).__name__ != "int" or type(row1end).__name__ != "int" or type(row2end).__name__ != "int":
    print(bcolors.FAIL+"   [ERROR] rows need to be integers "+bcolors.ENDC)
    print(bcolors.FAIL+"   [ERROR] ABORT"+bcolors.ENDC)
    sys.exit() 
    

print(bcolors.BOLD+"   Read Fits files" +bcolors.ENDC)
## open fits files and read the data sections into 2 dimensional arrays
# file with spectrum
spec1File     = fits.open('%s'%spec1FileName)
spec1FileData = spec1File[0].data
# darkframe for the spectrum
spec2File     = fits.open('%s'%spec2FileName)
spec2FileData = spec2File[0].data


## get dimensions of the data
dataColumns = len(spec1FileData[0])      # number of columns in the file
dataLines   = len(spec1FileData)         # number of lines in the file
N = dataColumns

## read lines that contain spectral information

print(bcolors.BOLD+"   Extract spectra " +bcolors.ENDC)
spec1 = list(np.average(spec1FileData[row1start:row1end],axis=0))
spec2 = list(np.average(spec2FileData[row2start:row2end],axis=0))
#spec1 = []
#tmp = 0
#for j in range(0, dataColumns):
    #for i in range(row1start,row1end):
        #tmp += spec1FileData[i][j]
    #tmp /= (row1end-row1start)
    #spec1.append(tmp)

#spec2 = []
#tmp = 0
#for j in range(0, dataColumns):
    #for i in range(row2start,row2end):
        #tmp = tmp + spec2FileData[i][j]
    #tmp = tmp / (row2end-row2start)
    #spec2.append(tmp)

#print(spec1,spec2)


#Set automatic plotrange
if len(spec1)+offset_spec1 > len(spec2)+offset_spec2:
    xmax = 1.02*(len(spec1)+offset_spec1)
else:
    xmax = 1.02*(len(spec2)+offset_spec2)

if offset_spec1 < offset_spec2:
    xmin = offset_spec1-.02*(len(spec1)+offset_spec2)
else:
    xmin = offset_spec2-.02*(len(spec1)+offset_spec1)


fig = plt.figure(figsize=(10,5))


print(bcolors.BOLD+"   Plot spectra "+bcolors.ENDC)


#adjust plot range if zoom is selected
if zoom == 'yes':
    try:
        float(x_start)
        float(x_end)
    except ValueError:
        print(bcolors.WARNING+"   [WARNING] Zoom requested but range is not a number"+bcolors.ENDC)
        print(bcolors.WARNING+"             -> automatic plot range used"+bcolors.ENDC)
    else:
        if x_start > x_end:
            a = x_end
            x_end = x_start
            x_start = a
        xmin,xmax = float(x_start),float(x_end)
        print(bcolors.BOLD+"   Zoom in on range ("+bcolors.OKBLUE+str(x_start)+","+str(x_end)+")"+bcolors.ENDC)

plt.xlim([xmin,xmax])

#apply the normalization by dividing by max value amd plot spectrum
if normalize == "yes":
    if shift_spec1 > shift_spec2:
        pylab.ylim([0,1.1+shift_spec1])
    else:
        pylab.ylim([0,1.1+shift_spec2])
    plt.plot(np.arange(offset_spec1,len(spec1)+offset_spec1),spec1/max(spec1)+shift_spec1,'r-',linewidth=.5,label=os.path.split(spec1FileName)[-1][:-4])
    plt.plot(np.arange(offset_spec2,len(spec2)+offset_spec2),spec2/max(spec2)+shift_spec2,'b-',linewidth=.5,label=os.path.split(spec2FileName)[-1][:-4])
    spec1 = spec1/max(spec1)
    spec2 = spec2/max(spec2)
else:
    plt.plot(np.arange(offset_spec1,len(spec1)+offset_spec1),spec1+shift_spec1,'r-',linewidth=.5)
    plt.plot(np.arange(offset_spec2,len(spec2)+offset_spec2),spec2+shift_spec2,'b-',linewidth=.5)

#Set output filenames
if filename == "" or filename == "?":
    print(bcolors.WARNING+"   [WARNING] No filename given, use default 'spectrum'"+bcolors.ENDC)
    filename = "spectrum"
    
if filetype != "pdf" and filetype != "png" and filetype != "eps" and filetype != "ps" and filetype != "svg":
    print(bcolors.WARNING+"   [WARNING] No or unknown filetype given, use default 'pdf'"+bcolors.ENDC)
    filetype = "pdf"


    
#plot legende
plt.legend(loc=3,ncol=1)
#adjust ticks to appear all around
plt.tick_params(top=True,right=True,which='both',direction='in')
plt.minorticks_on()




#os.system('rm tmpdata')
#set plot labels
plt.xlabel('Pixel')
if normalize == 'yes':
    plt.ylabel('Normalized intensity')
else:
    plt.ylabel('Intensity')

print(bcolors.BOLD+"   Write "+bcolors.OKBLUE+filename + '.' + filetype+bcolors.ENDC)
#save figure
plt.savefig(filename + '.' + filetype,format=filetype,bbox_inches='tight')

#write spectral data into data set
if two_files == 'no':
    print(bcolors.BOLD+'   Write output/'+bcolors.OKBLUE+os.path.split(spec1FileName)[-1][:-3]+"dat" +bcolors.ENDC)
    os.system("mkdir -p output")
    os.system('touch output/'+os.path.split(spec1FileName)[-1][:-3]+"dat")
    f = open('output/'+os.path.split(spec1FileName)[-1][:-3]+"dat", 'w')
    for i in range(0,len(spec1)):
        j = i + offset_spec1
        f.write('%i\t%f\n'%(j,spec1[i]))
    f.close()

    print(bcolors.BOLD+'   Write output/'+bcolors.OKBLUE+os.path.split(spec2FileName)[-1][:-3]+"dat" +bcolors.ENDC)
    os.system('touch output/'+os.path.split(spec2FileName)[-1][:-3]+"dat")
    f = open('output/'+os.path.split(spec2FileName)[-1][:-3]+"dat", 'w')
    for i in range(0,len(spec2)):
        j = i + offset_spec2
        f.write('%i\t%f\n'%(j,spec2[i]))
    f.close()
    
else:
    print(bcolors.BOLD+'   Write output/'+bcolors.OKBLUE+output_file_name +bcolors.ENDC)
    os.system("mkdir -p output")
    os.system('touch output/'+output_file_name)
    f = open('output/'+output_file_name, 'w')
    for i in range(0,len(spec2)):
        j = i + offset_spec2
        f.write('{:>5}'.format(j)+' '+'{:>25}'.format(str(spec1[i]))+' '+'{:>25}'.format(str(spec2[i]))+'\n')
    f.close()
    
    
    
print(bcolors.OKGREEN+"   Done"+bcolors.ENDC)
