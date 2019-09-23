#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 22 18:50:42 2019

@author: samuelfielder

Here lies all the fuctions that didn't make it.
"""

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