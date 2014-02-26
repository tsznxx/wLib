#!/usr/bin/python
#Last-modified: 25 Feb 2014 02:52:07 PM

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
import scipy
#from wLib.wBed import IO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

class test(object):
    ''' Statistics test. '''
    def ttest(a,b, axis=0, equal_var= True):
        '''
        Student test.
        
        Parameters:
            a,b: array_like
                The arrays must have the same shape, except in the dimension corresponding to axis (the first, by default).
            axis: int, optional
                Axis can equal None (ravel array first), or an integer (the axis over which to operate on a and b).
            equal_var: bool, optional
                If True (default), perform a standard independent 2 sample test that assumes equal population variances [R248]. If False, perform Welch’s t-test, which does not assume equal population variance [R249].
        
        Returns:
            t: float or array
                The calculated t-statistic.
            prob : float or array
                The two-tailed p-value.
        '''
        return scipy.stats.ttest_ind(a,b,axis,equal_var)
    ttest=staticmethod(ttest)
    def ranksums(a,b):
        '''
        Compute the Wilcoxon rank-sum statistic for two samples.
        
        The Wilcoxon rank-sum test tests the null hypothesis that two sets of measurements are drawn from the same distribution. The alternative hypothesis is that values in one sample are more likely to be larger than the values in the other sample.

        This test should be used to compare two samples from continuous distributions. It does not handle ties between measurements in x and y. For tie-handling and an optional continuity correction see scipy.stats.mannwhitneyu.

        Parameters:
            x,y : array_like
                The data from the two samples
        
        Returns :   
            z-statistic : float
                The test statistic under the large-sample approximation that the rank sum statistic is normally distributed
            p-value : float
                The two-sided p-value of the test
        '''
        return scipy.stats.ranksums(a,b)
    ranksums=staticmethod(ranksums)
    def wilcoxon(x, y=None, zero_method='wilcox', correction=False)
        '''
        Calculate the Wilcoxon signed-rank test.
        
        The Wilcoxon signed-rank test tests the null hypothesis that two related paired samples come from the same distribution. In particular, it tests whether the distribution of the differences x - y is symmetric about zero. It is a non-parametric version of the paired T-test.
        
        Parameters :    
            x : array_like
                The first set of measurements.
            y : array_like, optional
                The second set of measurements. If y is not given, then the x array is considered to be the differences between the two sets of measurements.
            zero_method : string, {“pratt”, “wilcox”, “zsplit”}, optional
                “pratt”:
                    Pratt treatment: includes zero-differences in the ranking process (more conservative)
                “wilcox”:
                    Wilcox treatment: discards all zero-differences
                “zsplit”:
                    Zero rank split: just like Pratt, but spliting the zero rank between positive and negative ones
            correction : bool, optional
                If True, apply continuity correction by adjusting the Wilcoxon rank statistic by 0.5 towards the mean value when computing the z-statistic. Default is False.
        
        Returns :   
            T : float
                The sum of the ranks of the differences above or below zero, whichever is smaller.
            p-value : float
                The two-sided p-value for the test.
        
        Note:
            Because the normal approximation is used for the calculations, the samples used should be large. A typical rule is to require that n > 20.
        '''
        return scipy.stats.wilcoxon(x,y,zero_method,correction)
    wilcoxon=staticmethod(wilcoxon)

class Utils(object):
    ''' Methods frequently used in statistics. '''
    def percentile(sortedv, p, interpolation=False):
        '''
        Return the item at the percentile of a sorted list.
        
        Parameters:
            sortedv: list or array
                sorted list or array
            p: float
                percentile of the list. Range from 0 to 1.
            interpolation: bool
                If false (default), return the nearest item. Otherwise, return the interpolated value of the two nearest items.
        
        Returns:
            item: not sure
                the item nearest to the pth percentile.
        '''
        if interpolation:
            s = len(sortedv)*p
            idxl = int(s)
            pl = s - idxl
            if pl != 0:
                idxl -= 1
            return sortedv[idxl]*(1-pl) + sortedv[idxl+1]*pl
        return sortedv[round(percent*len(sortedv))-1]
    percentile=staticmethod(percentile)
    def fivenum(v):
        '''
        Tukey's five number summary (minimum, lower-hinge, median, upper-hinge, maximum).
        
        Parameters:
            v: list or array
                a list of numbers. 
        
        Returns:
            fv: tuple
                (minimum, lower-hinge, median, upper-hinge, maximum). Hinge is defined as 1.5 IQR from median.
        '''
        l = len(v)-1
        nv = numpy.array(v)
        nv.sort()
        md = numpy.median(nv)
        q1 = percentile(nv, 0.25)
        q3 = percentile(nv, 0.75)
        whisker = 1.5*(q3 - q1)
        md = numpy.median(nv)
        return (nv[0], md-whisker, md, md+whisker, nv[1])
    fivenum=staticmethod(fivenum)
    
    def rank(x):
        ''' return the rank of x.'''
        array=numpy.array(x)
        temp = array.argsort()
        ranks = numpy.arange(len(array))[temp.argsort()]
        return ranks
    rank=staticmethod(rank)
    
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
    tmean=staticmethod(tmean)

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

