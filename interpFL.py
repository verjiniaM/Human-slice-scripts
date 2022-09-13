# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 13:33:14 2016

"""
import numpy as np


def findLevels(A, level, mode='rising', boxWidth=0, rangeSubset=None):
    """Function to find level crossings in an 1d numpy array.  Based on the Igor
    function FindLevel. 
    Can find rising and/or falling crossings, control with the 'mode' paramter.
    Returns an numpy array of all crossings found and the number of crossings
    :param A: 1d numpy array
    :param level: floating point to search for in A
    :param mode: optional string: mode specfication. one of 'rising', 'falling' or 'both'
    :param boxWidth: optional int for local boxcar smoothing
    :param rangeSubset: optional list of ints to limit the search
    :returns: tuple, a numpy array of level crossings and the number of crossings
    """
    assert mode in ('rising', 'falling', 'both'), 'traceManip.findLevels: Unknown mode \'%s\'' % mode

    if boxWidth is not 0:
        A = np.convolve(A, np.array([1]*boxWidth)/float(boxWidth))

    crossings = np.diff(np.sign(A-level), axis=0)
    
    if mode is 'rising':
        rising_points = np.where(crossings > 0)
        return rising_points[0], len(rising_points[0])
    elif mode is 'falling':
        falling_points = np.where(crossings < 0)
        return falling_points[0], len(falling_points[0])
    else:
        all_crossing_points = np.where(np.abs(crossings) > 0)
        return all_crossing_points, len(all_crossing_points)

def findLevels1d(A, level, mode='rising', boxWidth=0):
    return findLevelsNd(A, level, mode=mode, axis=0, boxWidth=boxWidth)

def findLevelsNd(A, level, mode='rising', axis=0, boxWidth=0):
    """Function to find level crossings in an Nd numpy array. 
    Can find rising and/or falling crossings, control with the 'mode' paramter.
    Returns a binary array of level crossings, with true elements right AFTER a crossing.
    NOTE THAT THIS RETURNS DIFFERENT VALUES THAN findLevels().  if you want to get a list of
    locations where the crossings occurs, then use the following syntax:
    levels = findLevelsNd(array, level)
    level_crossings_locations = levels.nonzero()
    number_of_level_crossings = len(level_crossing_locations[0])
    Often, the crossings are noisy.  You can use np.diff() and findLevelsNd() again to help yourself out.
    :param A: 1d numpy array
    :param level: floating point to search for in A
    :param mode: optional string: mode specfication. one of 'rising', 'falling' or 'both'
    :param axis: optional integer, specifies dimension
    :param boxWidth: optional int for local boxcar smoothing
    :returns: binary array of level crossing locations
    """
    assert mode in ('rising', 'falling', 'both'), 'traceManip.findLevels: Unknown mode \'%s\'' % mode

    if boxWidth is not 0:
        A = nd.convolve1d(A, np.array([1]*boxWidth)/float(boxWidth), axis=axis)

    crossings = np.diff(np.sign(A-level), axis=axis)
    
    if mode is 'rising':
        return crossings>0
    elif mode is 'falling':
        return crossings<0
    else:
        return np.abs(crossings>0)
