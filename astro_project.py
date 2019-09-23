#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 19:44:17 2019

@author: samuelfielder
"""


import numpy as np
from astropy.convolution import convolve_fft, Box2DKernel
import os

from definitions import galaxy_table_maker
from definitions import name_splitter
from definitions import data_grabber
from definitions import morph_finder
from definitions import reproject_fits
from definitions import string_trimmer
from definitions import box_maker
from definitions import data_visualizer

# For Script Timing
from datetime import datetime
startTime = datetime.now()

import yaml
# Import Necessary Parameters here
with open("config.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.BaseLoader)
    
IMAGE_DIR = cfg['parameters']['IMAGE_DIR']
CONV_SIZE = int(cfg['parameters']['CONV_SIZE'])
BOX_SIZE = int(cfg['parameters']['BOX_SIZE'])
SHOW_PLOTS = cfg['parameters']['SHOW_PLOTS']

# Creating Nested List of Filenames
data = galaxy_table_maker(IMAGE_DIR)

# Splitting Galaxy names off of directory strings
galaxies = []
for i in range(len(data[0])):
    name = name_splitter(data[0][i])
    galaxies.append(name)


if os.path.isfile('morphology_list.txt'):
    print('File Detected')
    with open('morphology_list.txt', 'r') as file:
        morphology = file.readlines()
else:
    print('File Not Detected')
    morphology = morph_finder(galaxies)
    with open('morphology_list.txt', 'w+') as file:
        for item in morphology:
            file.write("%s\n" % item)

# Label List Creation
galaxy_labels = []
for i in range(len(morphology)):    
    galaxy_labels.append(string_trimmer(morphology[i]))

#%%
""" Here starts the Main Loop for Pre-Processing"""
index = 0
phot_1_header, phot_1_data = data_grabber(data[0][index])
phot_2_header, phot_2_data = data_grabber(data[1][index])
wt_1_header, wt_1_data = data_grabber(data[2][index])
wt_2_header, wt_2_data = data_grabber(data[3][index])

box_kernel = Box2DKernel(CONV_SIZE)

# Reproject 2 Data/Wt in Header 1 Coordinates
phot_2_rep_data, _ = reproject_fits(data[1][index], phot_1_header)
wt_2_rep_data, _ = reproject_fits(data[3][index], wt_1_header)

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



print("Time for Script to Complete: ", datetime.now() - startTime)

