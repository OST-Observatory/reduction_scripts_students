#!/bin/python3

# BARYCENTRIC VELOCITY CORRECTION CALCULATOR v1.0

##################################################################################
############################ Script Parameters ###################################
##################################################################################
####################   REPLACE EVERY "?" WITH ITS RESPECTIVE VALUE    ###################


##### Target coordinates #####
# enter the object's ICRS coordinates

right_ascension                     = '?h ?m ?s'         #as hour angle: hour | minute | second
declination                         = '+?deg ?m ?s'      #in deg | arcmin |arcsec


##### Observation time #####

observation_date                    = '????-??-?? ??:??:??' #enter the observation date formatted as 'YYYY-MM-DD HH:MM:SS'. Be careful, to enter the correct day if observations were made past midnight!
time_zone                           = ?                    #type +1 for the central european time zone or +2 when time was switched to daylight saving



##### Observatory #####

observatory_latitude                = ?                    #values in degrees, north is positive
observatory_longitude               = ?                    #values in degrees, east is positive
observatory_height_above_MSL        = ?                    #enter in meter
##################################################################################




############################ Do not change anything below this line ###################################
##################################################################################

####################   IMPORT OF PACKAGES    ####################
import sys
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u
from astropy.units.function import units


####################   DEFINE COLOR FORMATTING    ####################
class bcolors:                                              #definition of some useful colorcodes for a color-formatted output in the terminal
    HEADER          = '\033[95m'
    OKBLUE          = '\033[94m'
    OKGREEN         = '\033[92m'
    WARNING         = '\033[93m'
    FAIL            = '\033[91m'
    GRAY            = '\033[90m'
    ENDC            = '\033[0m'
    BOLD            = '\033[1m'
    UNDERLINE       = '\033[4m'


####################   CONSTRUCTION OF ASTROPY SKYCOORD INSTANCE    ####################
try:
    int(time_zone)
except:
    print(
        bcolors.FAIL + '[FATAL ERROR]' + bcolors.ENDC + '   The time zone you have entered looks odd and can not be interpreted with this program. Make sure you entered an integer value e.g. +1 for CET or +2 for CEST.'
    )

try:
    obs_time        = Time(observation_date, format='iso', scale='utc') - time_zone*u.hour
except:
    print(
        bcolors.FAIL + '[FATAL ERROR]' + bcolors.ENDC + '   The observation time could not be identified. You entered: ' + observation_date \
        + '\nPlease check the date and time you entered. Make sure to use a format like "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD HH:MM"'
    )
    sys.exit()

try:    #construct an EarthLocation instance with the parameters given
    obs_site = EarthLocation(
        lat         = observatory_latitude,
        lon         = observatory_longitude,
        height      = observatory_height_above_MSL
    )
except:
    print(
        bcolors.FAIL + '[FATAL ERROR]' + bcolors.ENDC + '   Could not parse a location on earth from the values you entered.' \
        + '\nPlease check the lat., lon. and height values you entered. Make sure you use degrees and meters as units.'
    )
    sys.exit()

try:    #construct a SkyCoord instance with the parameters given
    obs_coord = SkyCoord(
        ra          = right_ascension,
        dec         = declination,
        unit        = (u.hourangle, u.deg),
        frame       = 'icrs',
        obstime     = obs_time,
        location    = obs_site
    )
except:
    print(
        bcolors.FAIL + '[FATAL ERROR]' + bcolors.ENDC + '   Could not create the SkyCoordinate to compute the barycentric velocity correction.' \
        + '\nPlease check the right ascension and declination values you entered. You entered:\n\tright ascension:   ' + right_ascension \
        + '\n\tdeclination:   ' + declination + '\nMake sure the right ascension is formatted like "??h ??m ??s" and for the declination "??deg ??m ??s".'
    )
    sys.exit()

####################   PRINT THE OUTPUT    ####################

if time_zone == 0:
    zone            = ' UTC'
elif time_zone == 1:
    zone            = ' CET'
elif time_zone == 2:
    zone            = ' CEST'
else:
    zone            = ' UTC ' + str(time_zone)

print('\n' + bcolors.OKGREEN + '[SUCCESS] ' + bcolors.ENDC + 'The barycentric velocity correction could be calculated using the following parameters:')

print('Date and Time of observation: ' + bcolors.OKBLUE + observation_date + zone + bcolors.ENDC)

print('\nThe ICRS coordinates of the observed object are: \n\t right ascension:   ' + bcolors.OKBLUE \
    + str(obs_coord.ra.to(u.hourangle)) + bcolors.ENDC + '\n\t declination:   ' \
    + bcolors.OKBLUE + str(obs_coord.dec)+ bcolors.ENDC)

print('\nThe object therefore is located in the constellation ' + bcolors.OKBLUE + str(obs_coord.get_constellation()) + bcolors.ENDC)

print('During its observation the object stood ' + bcolors.OKBLUE + str('{0:.2f}'.format(obs_coord.transform_to(AltAz()).alt.value)) + \
     ' deg ' + bcolors.ENDC + 'above the horizon')

bvc = obs_coord.radial_velocity_correction(kind='barycentric')
print('\n======================================================')
print('This leads to a barycentric velocity correction of: ' + '\n\t ' + r'v_{bvc} = ' + bcolors.OKGREEN + str(bvc.to(u.km / u.s)) + bcolors.ENDC)
print('======================================================\n\n')
