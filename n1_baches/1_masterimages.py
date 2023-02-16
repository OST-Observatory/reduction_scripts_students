#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
    Reduce spectral data taken with the BACHES spectrograph for extraction
    and analysis by means of MIDAS
'''

############################################################################
####           Configuration: modify the file in this section           ####
############################################################################

###
#   Path to the directories with the images
#
#   Darks:
path_darks = '?'

#   Flat darks:
path_flat_darks = '?'

#   Flats:
path_flats = '?'

#   Thorium Argon exposures:
path_thars = '?'

#   Spectra:
path_spectra = '?'

#   Output directory for the reduced flats. The master files will be saved in
#   the current directory.
out_path = 'output'


###
#   Flip images? Possibilities: True and False
#
flip_bool = False
# flip_bool = True


###
#   Bin the images? Possibilities: True and False
#
bin_bool = False
# bin_bool = True
#   Binning factor
bin_val = 2


###
#   Trim images to remove non essential parts and thus simplify MIDAS handling
#   Possibilities: True and False; Default: True
#
trim_bool = True

#   Number of pixel to be removed from the start (''_s'') and end (''_e'')
#   of the image in X (''_x_'') and Y (''_y_'') directory
#   Typically, the default settings below do not need to be changed!
trim_x_s=400
trim_x_e=400
trim_y_s=300
trim_y_e=250

############################################################################
####                            Libraries                               ####
############################################################################

import sys

from pathlib import Path

from astropy.stats import mad_std

import numpy as np

import ccdproc as ccdp

from astropy.nddata import CCDData
import astropy.units as u

from ost import checks, terminal_output
import ost.reduce.aux as aux

import warnings
warnings.filterwarnings('ignore')

from astropy import log
log.setLevel('ERROR')

############################################################################
####                            Functions                               ####
############################################################################

def master_img(path, out_path, img_type, flip_bool=False, bin_bool=False,
                bin_val=2, trim_bool=True, method='average',
                subtract_dark=False, master_dark='master_dark.fit',
                divide_flat=False, master_flat='master_flat.fit',
                scale_func=None, trim_x_s=400, trim_x_e=400, trim_y_s=300,
                trim_y_e=250):
    '''
        Create master images

        Parameters
        ----------
        path                : `string`
            Path to the directory with the files

        out_path            : `string`
            Path to the directory to which the individual files should
            be written

        img_type            : `string`
            String that characterizes the image type.

        flip_bool           : `boolean`, optional
            If `True` the images will be flipped in X and Y direction
            Default is ``False``.

        bin_bool            : `boolean`, optional
            If `True` the images will be binned in X and Y direction
            Default is ``False``.

        bin_val             : `integer`, optional
            Value by which the image should be reduced in size.
            Default is ``2``.

        trim_bool           : `boolean`, optional
            If `True` the images will be trimmed in X and Y direction
            Default is ``True``.

        method              : `string`, optional
            Method used to average the images.
            Default is ``average``.

        subtract_dark       : `boolean`, optional
            If `True` a master dark will be subtracted for all input images.
            Default is ``False``.

        master_dark         : `string`, optional
            Name of the master dark file to be subtracted from all input
            images.
            Default is ``master_dark.fit``.

        divide_flat         : `boolean`, optional
            If `True` the spectra will be divided by a master flat.
            Default is ``False``.

        master_flat         : `string`, optional
            Name of the master flat file.
            Default is ``master_flat.fit``.

        scale_func          : `function` or `None`, optional
            Return value of the provided function will be used to scale the
            image before combining. If `None` no scaling is applied.
            Default is ``None``.

        trim_x_s             : `integer`
            Number of Pixel to be removed from the start of the image in
            X direction.

        trim_x_e             : `integer`
            Number of Pixel to be removed from the end of the image in
            X direction.

        trim_y_s             : `integer`
            Number of Pixel to be removed from the start of the image in
            Y direction.

        trim_y_e             : `integer`
            Number of Pixel to be removed from the end of the image in
            Y direction.
    '''
    terminal_output.print_terminal(
        img_type,
        indent=1,
        string="Reduce {} images",
        )

    #   Load images
    images = ccdp.ImageFileCollection(path)

    #   Flip, trim, bin?
    if flip_bool:
        images = aux.flip_img(images, out_path / img_type)
    if bin_bool:
        images = aux.bin_img(images, out_path / img_type, bin_val)
    if trim_bool:
        images = aux.trim_img(
            images,
            out_path / img_type,
            xs=trim_x_s,
            xe=trim_x_e,
            ys=trim_y_s,
            ye=trim_y_e,
            )

    #   Subtract dark
    if subtract_dark:
        terminal_output.print_terminal(
            img_type,
            indent=2,
            string="Subtract master dark",
            )

        #   Read master dark
        dark = CCDData.read(master_dark)

        for img, file_name in images.ccds(
            ccd_kwargs={'unit': 'adu'},
            return_fname=True,
            ):

            #   Subtract the dark current
            img = ccdp.subtract_dark(
                img,
                dark,
                exposure_time='exptime',
                exposure_unit=u.second,
                )

            #   Save the result
            #directory = img_type+'_dark_corrected'
            img_path = out_path / img_type / 'dark_corrected'
            checks.check_out(img_path)
            img.write(img_path / file_name, overwrite=True)

        #   Reload images
        images = ccdp.ImageFileCollection(img_path)

    #   Divide by flat
    if divide_flat:
        terminal_output.print_terminal(
            img_type,
            indent=2,
            string="Divide by master flat",
            )

        #   Read master flat
        flat = CCDData.read(master_flat)

        for img, file_name in images.ccds(
            ccd_kwargs={'unit': 'adu'},
            return_fname=True,
            ):

            #   Divide by flat
            img = ccdp.flat_correct(img, flat)

            #   Save the result
            #directory = img_type+'_flatfielded'
            img_path = out_path / img_type / 'flatfielded'
            checks.check_out(img_path)
            img.write(img_path / file_name, overwrite=True)

        #   Reload images
        images = ccdp.ImageFileCollection(img_path)

    #   Apply filter to the image collection
    #   -> This is necessary so that the path to the image directory is
    #      added to the file names.
    images = images.filter(SIMPLE=True)

    #   Stack images
    combined_img = ccdp.combine(
        images.files,
        method=method,
        scale=scale_func,
        sigma_clip=True,
        sigma_clip_low_thresh=5,
        sigma_clip_high_thresh=5,
        sigma_clip_func=np.ma.median,
        sigma_clip_dev_func=mad_std,
        mem_limit=15e9,
        unit='adu',
        )

    #   Save master image
    combined_img.write('master_'+img_type+'.fit', overwrite=True)


############################################################################
####                               Main                                 ####
############################################################################

if __name__ == '__main__':
    ###
    #   Check input and output directories
    #
    path_darks = checks.check_pathlib_Path(path_darks)
    path_flat_darks = checks.check_pathlib_Path(path_flat_darks)
    path_thars = checks.check_pathlib_Path(path_thars)
    path_flats = checks.check_pathlib_Path(path_flats)
    if path_spectra != '?':
        path_spectra = checks.check_pathlib_Path(path_spectra)
    checks.check_out(out_path)
    out_path = Path(out_path)


    ###
    #   Master dark
    #
    master_img(
        path_darks,
        out_path,
        flip_bool=flip_bool,
        bin_bool=bin_bool,
        bin_val=bin_val,
        trim_bool=trim_bool,
        trim_x_s=trim_x_s,
        trim_x_e=trim_x_e,
        trim_y_s=trim_y_s,
        trim_y_e=trim_y_e,
        img_type='dark',
        )

    ###
    #   Master flat dark
    #
    master_img(
        path_flat_darks,
        out_path,
        flip_bool=flip_bool,
        bin_bool=bin_bool,
        bin_val=bin_val,
        trim_bool=trim_bool,
        trim_x_s=trim_x_s,
        trim_x_e=trim_x_e,
        trim_y_s=trim_y_s,
        trim_y_e=trim_y_e,
        img_type='flat_dark',
        )


    ###
    #   Master Thorium Argon
    #
    master_img(
        path_thars,
        out_path,
        flip_bool=flip_bool,
        bin_bool=bin_bool,
        bin_val=bin_val,
        trim_bool=trim_bool,
        trim_x_s=trim_x_s,
        trim_x_e=trim_x_e,
        trim_y_s=trim_y_s,
        trim_y_e=trim_y_e,
        img_type='thar',
        )


    ###
    #   Master flat
    #
    master_img(
        path_flats,
        out_path,
        flip_bool=flip_bool,
        bin_bool=bin_bool,
        bin_val=bin_val,
        trim_bool=trim_bool,
        trim_x_s=trim_x_s,
        trim_x_e=trim_x_e,
        trim_y_s=trim_y_s,
        trim_y_e=trim_y_e,
        img_type='flat',
        subtract_dark=True,
        master_dark='master_flat_dark.fit',
        scale_func=aux.inv_median,
        )


    ###
    #   Master spectrum
    #
    master_img(
        path_spectra,
        out_path,
        flip_bool=flip_bool,
        bin_bool=bin_bool,
        bin_val=bin_val,
        trim_bool=trim_bool,
        trim_x_s=trim_x_s,
        trim_x_e=trim_x_e,
        trim_y_s=trim_y_s,
        trim_y_e=trim_y_e,
        method='median',
        img_type='spectrum',
        subtract_dark=True,
        master_dark='master_dark.fit',
        #divide_flat=True,
        #master_flat='master_flat.fit',
        )
