#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 19:44:17 2019

@author: samuelfielder
"""


import numpy as np
import matplotlib.pyplot as plt
import astropy.wcs as wcs
from astropy.convolution import convolve, convolve_fft, Box2DKernel

from definitions import galaxy_table_maker
from definitions import galaxy_query
from definitions import name_splitter
from definitions import data_grabber
from definitions import morph_finder
from definitions import max_rectangle
from definitions import reproject_fits
from definitions import image_visualizer
from definitions import data_preparer
from definitions import string_trimmer
from definitions import box_maker

import logging
# Turns off logging for wcs coordinate missing -SIP
logging.getLogger("wcs").setLevel(logging.WARNING)

from datetime import datetime
startTime = datetime.now()

# Setting Static Variables
IMAGE_DIR = '/home/samuelfielder/Desktop/Astro_Project/Images/'
BOX_SIZE = 400

# Creating Nested List of Filenames
data = galaxy_table_maker(IMAGE_DIR)

# Splitting Galaxy names off of directory strings
galaxies = []
for i in range(len(data[0])):
    name = name_splitter(data[0][i])
    galaxies.append(name)


"""
# Finding Morphology for Galaxies
morphology = morph_finder(galaxies)
"""

# Temporary Import of Morphologies
with open('morphology_list.txt', 'r') as file:
    lines = file.readlines()

# Label List Creation
galaxy_labels = []
for i in range(len(lines)):    
    galaxy_labels.append(string_trimmer(lines[i]))
    

    
#test_data_1, test_header_1, test_data_2, \
# test_header_2 = data_preparer(data[0][0], data[1][0])
 
_, test_header_1, test_data_1 = data_grabber(data[0][2])

_, test_header_1_wt, test_data_1_wt = data_grabber(data[2][2])
 

box_kernel = Box2DKernel(275)
kernel_array = convolve_fft(test_data_1_wt, box_kernel, nan_treatment='fill')

# Returns (y,x) value of the array
ind = np.unravel_index(np.argmax(kernel_array, axis=None),
                       kernel_array.shape)

x_min, x_max, y_min, y_max = box_maker(ind, BOX_SIZE)
#%%
test_data_1_wt_resize = test_data_1_wt[y_min:y_max, x_min:x_max]

#%%
#image_visualizer(test_data_1, test_header_1)
#image_visualizer(test_data_1_wt, test_header_1_wt)
plt.imshow(test_data_1_wt)
plt.plot(ind[1], ind[0], 'r.')
plt.plot(x_min, y_min, 'r.')
plt.plot(x_min, y_max, 'r.')
plt.plot(x_max, y_min, 'r.')
plt.plot(x_max, y_max, 'r.')
plt.show()
plt.imshow(kernel_array)
plt.plot(ind[1], ind[0], 'r.')
#image_visualizer(test_data_2, test_header_2)


    
    
print("Time for Script to Complete: ", datetime.now() - startTime)

