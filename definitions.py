#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:13:07 2019

@author: samuelfielder

Main functions used in Astro_CNN.
"""

from astroquery.simbad import Simbad
from astropy.io import fits
import glob
from reproject import reproject_interp
import matplotlib.pyplot as plt
import astropy.wcs as wcs
import yaml

# To Turn off Warning Messages to console, from wcs library
import logging
# Turns off logging for wcs coordinate missing -SIP
# Logger name for WCS package is astropy
logging.getLogger("astropy").setLevel(logging.WARNING)

# Import Necessary Parameters here
with open("config.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile, Loader=yaml.BaseLoader)


def galaxy_dict_maker(directory):
    """Takes input directory and creates a dictionary with image data.

    Parameters
    ----------
    directory : string
        full directory string of fits file

    Returns
    -------
    data_dict : dictionary
        Contains phot_1, phot_2, wt_1, wt_2, mask_1, mask_2 objects.

    Notes:
        - Need to import config.ylm for necessary dictionary.
    """
    phot_1 = glob.glob(directory+cfg['search_strings']['phot_1'])
    phot_1.sort()

    phot_2 = glob.glob(directory+cfg['search_strings']['phot_2'])
    phot_2.sort()

    phot_1_wt = glob.glob(directory+cfg['search_strings']['wt_1'])
    phot_1_wt.sort()

    phot_2_wt = glob.glob(directory+cfg['search_strings']['wt_2'])
    phot_2_wt.sort()

    final_mask_1 = glob.glob(directory+cfg['search_strings']['mask_1'])
    final_mask_1.sort()

    final_mask_2 = glob.glob(directory+cfg['search_strings']['mask_2'])
    final_mask_2.sort()

    data_dict = {}

    data_dict['phot_1'] = phot_1
    data_dict['phot_2'] = phot_2
    data_dict['1_wt'] = phot_1_wt
    data_dict['2_wt'] = phot_2_wt
    data_dict['1_mask'] = final_mask_1
    data_dict['2_mask'] = final_mask_2

    return data_dict


def galaxy_query(name):
    """Queries astropy for the morphology of the given galaxy.

    Parameters
    ----------
    name : string
        galaxy name

    Returns
    -------
    morphology : string

    """
    custom = Simbad()
    custom.add_votable_fields('morphtype')
    result = custom.query_object(name)
    morphology = result[0]['MORPH_TYPE']

    return morphology


def name_splitter(name):
    """Takes string name of director, and splits off the galaxy name.

    Parameters
    ----------
    name : string

    Returns
    -------
    galaxy : string

    """
    galaxy = name.split('/')[-1].split('.')[0]

    return galaxy


def data_grabber(directory_string):
    """Takes a file location and extracts the image data alongside
        header and info on file.

    Parameters
    ----------
    directory_string : string

    Returns
    -------
    header : astropy.io.fits.header.Header
    data : numpy.ndarray

    """
    with fits.open(directory_string) as hdul:
        data = hdul[0].data
        header = hdul[0].header

    return header, data


def morph_finder(galaxy_list):
    """Finds morphological type of galaxy based on Astropy Query, with an
        galaxy name as an input.

    Parameters
    ----------
    galaxy_list : list
        list of strings
    Returns
    -------
    morph_list : list
        list of strings
    """
    morph_list = []
    for i in range(len(galaxy_list)):
        morph = galaxy_query(galaxy_list[i])
        morph_list.append(morph)

    return morph_list


def location_query(name):
    """Queries Astropy Simbad library for RA,DEC of given galaxy.

    Parameters
    ----------
    name : string

    Returns
    -------
    ra : string
    dec : string

    """
    results_table = Simbad.query_object(name)
    ra = results_table['RA']
    dec = results_table['DEC']

    return ra, dec


def reproject_fits(filename, header):
    """Reproject one data array from fits, into the WCS coordinates
        of another header.

    Parameters
    ----------
    filename : string
        Points to a fits file.
    header : astropy.io.fits.header.Header
        Header that filename data will be formatted to.

    Returns
    -------
    data : numpy.ndarray
    footprint : numpy.ndarray

    """
    with fits.open(filename) as hdu:
        data, footprint = reproject_interp(hdu[0], header)

    return data, footprint


def data_visualizer(phot_1_data, phot_1_header,
                    wt_1_data, wt_1_header,
                    conv_1,
                    phot_2_rep_data,
                    wt_2_rep_data,
                    conv_2_rep,
                    phot_2_data, phot_2_header,
                    wt_2_data, wt_2_header,
                    conv_2,
                    ind_phot_1=None,
                    ind_phot_2_rep=None):
    """Visualizes all important graphs in a 3x3 configuration.

    Parameters
    ----------
    phot_1_data : numpy.ndarray
    phot_1_header : astropy.io.fits.header.Header
    wt_1_data : numpy.ndarray
    wt_1_header : astropy.io.fits.header.Header
    conv_1 : numpy.ndarray
    phot_2_rep_data : numpy.ndarray
    wt_2_rep_data : numpy.ndarray
    conv_2_rep : numpy.ndarray
    phot_2_data : numpy.ndarray
    phot_2_header : astropy.io.fits.header.Header
    wt_2_data : numpy.ndarray
    wt_2_header : astropy.io.fits.header.Header
    conv_2 : numpy.ndarray
    ind_phot_1 : tuple
        Coordinates of np.arg(phot_1_data).
    ind_phot_2_rep : tuple
        Coordinates of np.arg(phot_2_rep_data).

    Returns
    -------

    Notes:
        - Will not return anything explicitly, but does invoke the command
        plt.show(), which will print graph to console.
    """

    vir = plt.cm.viridis

    fig = plt.figure(figsize=(10, 10))

    # Phot_1_data
    ax1 = fig.add_subplot(3, 3, 1, projection=wcs.WCS(phot_1_header))
    ax1.imshow(phot_1_data, cmap=vir)
    ax1.set_title('Phot 1')

    # wt_1 data
    ax2 = fig.add_subplot(3, 3, 2, projection=wcs.WCS(wt_1_header))
    ax2.imshow(wt_1_data, cmap=vir)
    ax2.set_title('Phot 1 wt')

    # conv_1 data
    ax3 = fig.add_subplot(3, 3, 3, projection=wcs.WCS(wt_1_header))
    ax3.imshow(conv_1)
    if ind_phot_1 is not None:
        ax3.plot(ind_phot_1[1], ind_phot_1[0], 'r.')
    ax3.set_title('Conv Phot 1 wt')

    # phot_2_rep data
    ax4 = fig.add_subplot(3, 3, 4, projection=wcs.WCS(phot_1_header))
    ax4.imshow(phot_2_rep_data, cmap=vir)
    ax4.set_title('Phot 2 Rep')

    # wt_2_rep data
    ax5 = fig.add_subplot(3, 3, 5, projection=wcs.WCS(wt_1_header))
    ax5.imshow(wt_2_rep_data, cmap=vir)
    ax5.set_title('Phot 2 Rep wt')

    # conv_2_rep data
    ax6 = fig.add_subplot(3, 3, 6, projection=wcs.WCS(wt_1_header))
    ax6.imshow(conv_2_rep)
    if ind_phot_2_rep is not None:
        ax6.plot(ind_phot_2_rep[1], ind_phot_2_rep[0], 'r.')
    ax6.set_title('Conv Phot 2 Rep wt')

    # phot_2 data
    ax7 = fig.add_subplot(3, 3, 7, projection=wcs.WCS(phot_2_header))
    ax7.imshow(phot_2_data, cmap=vir)
    ax7.set_title('Phot 2')

    # wt_2 data
    ax8 = fig.add_subplot(3, 3, 8, projection=wcs.WCS(wt_2_header))
    ax8.imshow(wt_2_data, cmap=vir)
    ax8.set_title('Phot 2 wt')

    # conv_2 data
    ax9 = fig.add_subplot(3, 3, 9, projection=wcs.WCS(wt_2_header))
    ax9.imshow(conv_2)
    ax9.set_title('Conv Phot 2 wt')

    # Show Plot
    plt.show()

    return


def galaxy_visualizer(phot_1, phot_2, vmin=0, vmax=1):
    """Plots

    Parameters
    ----------
    phot_1 : numpy.ndarray
    phot_2 : numpy.ndarray
    vmin : float
        Default : 0
        This is the lower limit for cmap argument.
    vmax : float
        Default : 1
        This is the higher limit for cmap argument.

    Returns
    -------

    Notes:
        - Will not return anything explicitly, but does invoke the command
        plt.show(), which will print graph to console.

    """
    # Setting the Colormap typing
    colormap = plt.cm.viridis
    fig = plt.figure(figsize=(10, 10))

    # Phot 1 Data
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.imshow(phot_1, cmap=colormap, vmin=vmin, vmax=vmax)
    ax1.set_title('Channel 1 Data')

    # Phot 2 Data
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.imshow(phot_2, cmap=colormap, vmin=vmin, vmax=vmax)
    ax2.set_title('Channel 2 Data')

    plt.show()

    return


def string_trimmer(morph):
    """Returns the single typing of the morphological input string.

    Parameters
    ----------
    morph : string
        Full string returned by astropy.Simbad query.

    Returns
    -------
    morph_type : string
        Single letter representing chosen morphological type.

    Notes:
        - This function will detect: 'S', 'E', and 'S0' morph typing,
        and set them according to a pre defined layout.

    """

    morph_type = str()

    if 'E' in morph:
        morph_type = str('E')
    elif 'S' in morph:
        if 'S0' in morph:
            morph_type = str('E')
        else:
            morph_type = str('S')
    else:
        morph_type = str(None)

    return morph_type


def box_maker(center_indexes, BOX_SIZE):
    """Returns new bounds for box, around center index given.

    Parameters
    ----------
    center_indexes : tuple
        Due to convention, x and y are flipped arguments.
    BOX_SIZE : int
        This is the length of one side of the box.

    Returns
    -------
    x_min : int
    x_max : int
    ymin : int
    y_max : int

    """

    x_min = int(center_indexes[1] - (BOX_SIZE/2))
    x_max = int(center_indexes[1] + (BOX_SIZE/2))
    y_min = int(center_indexes[0] - (BOX_SIZE/2))
    y_max = int(center_indexes[0] + (BOX_SIZE/2))

    return x_min, x_max, y_min, y_max
