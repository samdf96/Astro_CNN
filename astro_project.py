#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 19:44:17 2019

@author: samuelfielder
"""

import numpy as np
from astropy.convolution import convolve_fft, Box2DKernel
import os
from astropy.io import fits
import glob
import yaml

from definitions import galaxy_dict_maker
from definitions import name_splitter
from definitions import data_grabber
from definitions import morph_finder
from definitions import reproject_fits
from definitions import string_trimmer
from definitions import box_maker
from definitions import data_visualizer
from definitions import galaxy_visualizer

# For Script Timing
from datetime import datetime
startTime = datetime.now()

# Import Necessary Parameters here
with open("config.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.BaseLoader)

IMAGE_DIR = cfg['parameters']['IMAGE_DIR']
SAVE_DIR = cfg['parameters']['SAVE_DIR']

# Creating SAVE_DIR and subdirectories if none found
if not os.path.isdir(SAVE_DIR):
    os.mkdir(SAVE_DIR)
if not os.path.isdir(SAVE_DIR + 'current/'):
    os.mkdir(SAVE_DIR + 'current/')
if not os.path.isdir(SAVE_DIR + 'backup/'):
    os.mkdir(SAVE_DIR + 'backup/')

# Creating Current and Backup locations as callable strings
CURRENT_DIR = SAVE_DIR + 'current/'
BACKUP_DIR = SAVE_DIR + 'backup/'

# Setting Parameters
CONV_SIZE = int(cfg['parameters']['CONV_SIZE'])
BOX_SIZE = int(cfg['parameters']['BOX_SIZE'])
# Converting from str() to bool()
if cfg['parameters']['SHOW_PLOTS'] == 'False':
    SHOW_PLOTS = False
else:
    SHOW_PLOTS = True

# Setting Skips for Galaxy Loop
skips = cfg['parameters']['GAL_SKIP']
skips = [int(i) for i in skips]

# Creating Nested List of Filenames
data = galaxy_dict_maker(IMAGE_DIR)

# Splitting Galaxy names off of directory strings
galaxies = []
for i in range(len(data['phot_1'])):
    name = name_splitter(data['phot_1'][i])
    galaxies.append(name)


if os.path.isfile('morphology_list.txt'):
    print('Morphology_list.txt file detected. Reading in data.')
    with open('morphology_list.txt', 'r') as file:
        morphology = file.readlines()
else:
    print('Morphology_list.txt file not detected. Querying Simbad,\
     and writing to file.')
    morphology = morph_finder(galaxies)
    with open('morphology_list.txt', 'w+') as file:
        for item in morphology:
            file.write("%s\n" % item)


# Setting empty list for final data output
phot_data = []
label_data = []

""" Here starts the Main Loop for Pre-Processing"""
for index in range(len(data['phot_1'])):
    print('Working on index: ', index, '/', len(data['phot_1']))
    
    # Skip Galaxies that are not working correctly
    if index in skips:
        continue
    phot_1_header, phot_1_data = data_grabber(data['phot_1'][index])
    phot_2_header, phot_2_data = data_grabber(data['phot_2'][index])
    wt_1_header, wt_1_data = data_grabber(data['1_wt'][index])
    wt_2_header, wt_2_data = data_grabber(data['2_wt'][index])

    box_kernel = Box2DKernel(CONV_SIZE)

    # Reproject 2 Data/Wt in Header 1 Coordinates
    phot_2_rep_data, _ = reproject_fits(data['phot_2'][index], phot_1_header)
    wt_2_rep_data, _ = reproject_fits(data['2_wt'][index], wt_1_header)

    # Convolution for all wt images
    conv_1 = convolve_fft(wt_1_data,
                          box_kernel,
                          nan_treatment='fill')
    conv_2 = convolve_fft(wt_2_data,
                          box_kernel,
                          nan_treatment='fill')
    conv_2_rep = convolve_fft(wt_2_rep_data,
                              box_kernel,
                              nan_treatment='fill')

    # Finding Max Coordinate in phot_1 and phot_2_rep images
    ind_phot_1 = np.unravel_index(np.argmax(conv_1, axis=None),
                                  conv_1.shape)

    ind_phot_2_rep = np.unravel_index(np.argmax(conv_2_rep, axis=None),
                                      conv_2_rep.shape)

    if SHOW_PLOTS:
        data_visualizer(phot_1_data, phot_1_header,
                        wt_1_data, wt_1_header, conv_1,
                        phot_2_rep_data, wt_2_rep_data, conv_2_rep,
                        phot_2_data, phot_2_header,
                        wt_2_data, wt_2_header, conv_2,
                        ind_phot_1,
                        ind_phot_2_rep)

    # Creating actual Phot_1 Dataset
    xmin1, xmax1, ymin1, ymax1 = box_maker(ind_phot_1, BOX_SIZE)
    phot_1_data_res = phot_1_data[ymin1:ymax1, xmin1:xmax1]

    xmin2, xmax2, ymin2, ymax2 = box_maker(ind_phot_2_rep, BOX_SIZE)
    phot_2_rep_data_res = phot_2_rep_data[ymin2:ymax2, xmin2:xmax2]

    if SHOW_PLOTS:
        galaxy_visualizer(phot_1_data_res, phot_2_rep_data_res)

    # Flattening Arrays and combining into one large array
    phot_1_flat = np.array(phot_1_data_res).flatten()
    phot_2_flat = np.array(phot_2_rep_data_res).flatten()

    phot_data.append(np.concatenate((phot_1_flat, phot_2_flat)))
    label_data.append(string_trimmer(morphology[index]))


# FITS File creation and data saving

format_string_data = str(len(phot_data[0])) + 'D'

data_column = fits.Column(name='data',
                          format=format_string_data,
                          array=phot_data)
label_column = fits.Column(name='label',
                           format='20A',
                           array=label_data)
galaxy_column = fits.Column(name='galaxy',
                            format='20A',
                            array=galaxies)
morphology_column = fits.Column(name='morphology',
                                format='20A',
                                array=morphology)

# Creating ColDefs object
cols = fits.ColDefs([galaxy_column,
                     morphology_column,
                     data_column,
                     label_column])

hdu = fits.BinTableHDU.from_columns(cols)

if len(glob.glob(CURRENT_DIR + '*.fits')) > 0:
    file_name = glob.glob(CURRENT_DIR + '*.fits')[0]
    fits_name = file_name.split('/')[-1]
    os.rename(file_name, BACKUP_DIR + fits_name)

save_time = f'{datetime.now():%Y-%m-%d_%H:%M:%S}'
fits_save = CURRENT_DIR + 'processed_data_' + save_time + '.fits'
hdu.writeto(fits_save)

print("Time for Script to Complete: ", datetime.now() - startTime)
