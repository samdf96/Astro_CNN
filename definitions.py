#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:13:07 2019

@author: samuelfielder

Main functions used in project.
"""

from astroquery.simbad import Simbad
from astropy.io import fits
import glob
from reproject import reproject_interp
import matplotlib.pyplot as plt
import astropy.wcs as wcs


def galaxy_table_maker(directory):
    """ Takes input directory and creates table with six columns, each with
        the following: Ch1, Ch2, Ch1_wt, Ch2_wt, Ch1_mask, Ch2_mask.

        Inputs:
            - directory: string
        Outputs:
            - image_lists: nested list
    """
    phot_1 = glob.glob(directory+'*.phot.1.fits')
    phot_1.sort()

    phot_2 = glob.glob(directory+'*.phot.2.fits')
    phot_2.sort()

    phot_1_wt = glob.glob(directory+'*.phot.1_wt.fits')
    phot_1_wt.sort()

    phot_2_wt = glob.glob(directory+'*.phot.2_wt.fits')
    phot_2_wt.sort()

    final_mask_1 = glob.glob(directory+'*.1.final_mask.fits')
    final_mask_1.sort()

    final_mask_2 = glob.glob(directory+'*.2.final_mask.fits')
    final_mask_2.sort()

    image_lists = []

    image_lists.append(phot_1)
    image_lists.append(phot_2)
    image_lists.append(phot_1_wt)
    image_lists.append(phot_2_wt)
    image_lists.append(final_mask_1)
    image_lists.append(final_mask_2)

    return image_lists


def galaxy_query(name):
    """ Querys astropy for the morphology of the given galaxy.

    Inputs:
        - name: string
    Outputs:
        - morphology: string
    """
    custom = Simbad()
    custom.add_votable_fields('morphtype')
    result = custom.query_object(name)
    morphology = result[0]['MORPH_TYPE']

    return morphology


def name_splitter(name):
    """ Takes string name of directory, and splits off the galaxy name.

    Inputs:
        - name: string
    Output:
        - galaxy: string
    """
    galaxy = name.split('/')[-1].split('.')[0]

    return galaxy


def data_grabber(directory_string):
    """ Takes a file location and extracts the image data alongside
        header and info on file.

    Inputs:
        - directory_string: string
    Outputs:
        - info: bound method (contains text)
        - header: astropy.io.fits.header.Header
        - data: numpy.ndarray
    """
    with fits.open(directory_string) as hdul:
        info = hdul.info
        data = hdul[0].data
        header = hdul[0].header

    return info, header, data


def morph_finder(galaxy_list):
    """ Finds morphological type of galaxy based on Astropy Query, with an
        galaxy name as an input.

    Inputs:
        - galaxy_list: list of strings
    Outputs:
        - morph_list: list of strings
    """
    morph_list = []
    for i in range(len(galaxy_list)):
        morph = galaxy_query(galaxy_list[i])
        morph_list.append(morph)

    return morph_list


def location_query(name):
    results_table = Simbad.query_object(name)
    ra = results_table['RA']
    dec = results_table['DEC']

    return ra, dec


def kadanes(array):
    """ Finds sum of maximum subarray of inputted array. Also returns
        the starting and ending indexes of the found subarray.

    Inputs:
        - array: numpy.array
    Outputs:
        - maximum: float
        - max_start: float
        - max_end: float
    """
    # Initialize Values
    maximum = 0
    max_start = -1
    max_end = -1
    current_start = 0
    ongoing_max = 0

    for i in range(0, len(array)):
        ongoing_max += array[i]

        if ongoing_max < 0:
            ongoing_max = 0
            current_start = i + 1

        if ongoing_max > maximum:
            max_start = current_start
            max_end = i
            maximum = ongoing_max

    return maximum, max_start, max_end


def max_rectangle(array):
    """ Finds the maximum sum subarray of any inputted array.
        Uses the kadanes definition for sub-routines.

    Inputs:
        - array: numpy.array
    Outputs:
        - rec_max: float
            - Sum of found subarray
        - (rec_left, rec_right): tuple
            - Right to left indexes of found subarray
        - (rec_top, rec_bot): tuple
            - Top to bottom indexes of found subarray
    """
    # Setting Parameters
    rows = len(array)
    cols = len(array[0])

    # Initialize Values
    max_sum = float("-inf")
    rec_left, rec_right, rec_top, rec_bot = (-1, -1, -1, -1)

    # Looping through all iterations of columns
    for left in range(cols):
        temp = [0 for _ in range(rows)]  # Setting temp array

        for right in range(left, cols):

            # Iterating through all rows
            for i in range(rows):
                temp[i] += array[i][right]  # Summing for Kadanes

            maximum, max_start, max_end = kadanes(temp)
            if maximum > max_sum:
                max_sum = maximum
                rec_left = left
                rec_right = right
                rec_top = max_start
                rec_bot = max_end

    return max_sum, (rec_left, rec_right), (rec_top, rec_bot)


def reproject_fits(filename, header):
    """ Reproject one data array from fits, into the WCS coordinates
        of another header. This will return a footprint array, as well
        as a data array.

    Inputs:
        - hdu:
        - header: astropy.io.fits.header.Header
    Ouput:
        - data: numpy.array
        - footprint: numpy.array
    """
    with fits.open(filename) as hdu:
        data, footprint = reproject_interp(hdu[0], header)

    return data, footprint


def image_visualizer(data, header, plot_show=True, plot_return=False):
    """ Matplotlib plot using WCS coordinates. User can choose to
        return (fig,ax) components, and/or for the plot to be shown
        automatically.

    Inputs:
        - data: numpy.array
        - header: astropy.io.fits.header.Header
        - plot_show: boolean
        - plot_return: boolean
    """
    # Creation of the WCS object
    wcs_obj = wcs.WCS(header)
    fig = plt.figure()
    fig.add_subplot(111, projection=wcs_obj)
    plt.imshow(data, origin='lower', cmap=plt.cm.viridis)
    plt.xlabel('RA')
    plt.ylabel('Dec')

    if plot_show:
        plt.show()

    if plot_return:
        return fig

    return


def data_preparer(filename_ch1, filename_ch2):
    """ Wrapper function that calls both channel images and returns
        both data arrays, with the Ch2 reprojected into Ch1 Header
        coordinates.

    Inputs:
        - filename_ch1: string
        - filename_ch2: string
    Outputs:
        - data_ch1: numpy.array
        - data_ch2_rep: numpy.array
    """
    _, header_ch1, data_ch1 = data_grabber(filename_ch1)
    _, header_ch2, data_ch2 = data_grabber(filename_ch2)

    # Reprojection Ch2 in terms of Ch1 Header

    data_ch2_rep, _ = reproject_fits(filename_ch2, header_ch1)

    return data_ch1, header_ch1, data_ch2_rep, header_ch2


def string_trimmer(morph):
    """ Returns singular typing of morphology."""
    # Creating list of characters from input
    morph_type = str()
    morph_char = list(morph)

    if ("E" in morph_char) & ("S" not in morph_char):
        morph_type == str("E")
    elif ("S" in morph_char) & ("E" not in morph_char):
        morph_type = str("S")
    elif ("S" in morph_char) & ("E" in morph_char):
        morph_type = str("E")
    else:
        morph_type == None

    return morph_type


def box_maker(center_indexes, BOX_SIZE):
    """ Returns new bounds for box, around center index."""
    
    x_min = int(center_indexes[1] - (BOX_SIZE/2))
    x_max = int(center_indexes[1] + (BOX_SIZE/2))
    y_min = int(center_indexes[1] - (BOX_SIZE/2))
    y_max = int(center_indexes[1] + (BOX_SIZE/2))
    
    return x_min, x_max, y_min, y_max
