#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 19:44:17 2019

@author: samuelfielder
"""


import numpy as np
import matplotlib.pyplot as plt
import astropy.wcs as wcs

from definitions import galaxy_table_maker
from definitions import galaxy_query
from definitions import name_splitter
from definitions import data_grabber
from definitions import morph_finder
from definitions import max_rectangle
from definitions import reproject_fits
from definitions import image_visualizer
from definitions import data_preparer

import logging
# Turns off logging for wcs coordinate missing -SIP
logging.getLogger("wcs").setLevel(logging.WARNING)

IMAGE_DIR = '/home/samuelfielder/Desktop/Astro_Project/Images/'

#Creating Nested List of Filenames
data = galaxy_table_maker(IMAGE_DIR)

#Splitting Galaxy names off of directory strings
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

#%%

test_data_1, test_header_1, test_data_2, test_header_2 = data_preparer(data[0][0], data[1][0])


#%%
image_visualizer(test_data_1, test_header_1)
image_visualizer(test_data_2, test_header_2)
#%%
#test_wt_data = np.nan_to_num(test_wt_data, nan=1)
#test_wt_data[np.isnan(test_wt_data)]=float("-inf")

#start = time.time()
#result = max_rectangle(test_wt_data)
#print("Elapsed time for Max Rectangle: ", time.time() - start)