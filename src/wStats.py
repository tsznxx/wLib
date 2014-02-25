#!/usr/bin/python
#Last-modified: 19 Oct 2013 03:37:27 PM

#         Module/Scripts Description
# 
# Copyright (c) 2008 Yunfei Wang <Yunfei.Wang1@utdallas.edu>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: Yunfei.Wang1@utdallas.edu

# ------------------------------------
# python modules
# ------------------------------------

import os,sys
import numpy
import string
#from wLib.wBed import IO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def percentile(sortedv, percent):
    ''' return the percentile of a sorted list. '''
    return sortedv[round(percent*len(sortedv))-1]

def fivenum(v):
    """Returns Tukey's five number summary (minimum, lower-hinge, median, upper-hinge, maximum) for the input vector, a list or array of numbers based on 1.5 times the interquartile distance"""
    l = len(v)-1
    nv = numpy.array(v)
    nv.sort()
    md = numpy.median(nv)
    q1 = percentile(nv, 0.25)
    q3 = percentile(nv, 0.75)
    whisker = 1.5*(q3 - q1)
    md = numpy.median(nv)
    return (nv[0], md-whisker, md, md+whisker, nv[1])

def rank(x):
    ''' return the rank of x.'''
    array=numpy.array(x)
    temp = array.argsort()
    ranks = numpy.arange(len(array))[temp.argsort()]
    return ranks

def tmean(x,lower=0,upper=1):
    ''' 
        trimmedmean of x. 
        Example: tmean (x,0.05,0.95)    
    '''
    array=numpy.array(x)
    if lower==0 and upper==1:
        return numpy.mean(array)
    else:
        array.sort()
        return numpy.mean(array[int(array.size*lower):int(array.size*upper)+1])

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

