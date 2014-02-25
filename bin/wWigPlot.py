#!/usr/bin/python
#Last-modified: 27 Oct 2013 10:36:31 PM

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

from pylab import *
from wLib import IO

# ------------------------------------
# constants
# ------------------------------------

cdict  = {'red':  ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 1.0, 1.0),
                   (1.0, 0.8, 0.8)),

         'green': ((0.0, 0.0, 0.0),
                   (0.1, 0.0, 0.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.8, 0.8),
                   (0.1, 1.0, 1.0),
                   (0.5, 1.0, 1.0),
                   (0.9, 0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }

# ------------------------------------
# Misc functions
# ------------------------------------
def smooth(x,window_len=11,window='hanning'):
    '''Input is a matrix. Smooth the depth for each interval.'''
    nrow, ncol =x.shape
    if ncol < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
    #s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=mp.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    w/=w.sum()
    y = np.apply_along_axis( lambda row:np.convolve(w,row,mode='same'),1,x)
    return y

def standarize(data):
    # Data processing
    # remove outlier
    data = np.minimum(data, np.mean(data)+3*np.std(data)) # data > mean+3*std will be mean+3*std
    # smooth data
    data = smooth(data,11)
    # standarize
    data = (data-np.mean(data))/(np.max(data)-np.min(data))
    return data
    
def Plot(data, xticks=[], yticks=[], order=[], showColBar=False):
    nrow,ncol = data.shape
    # change the order of data
    if len(order) >0:
        data = data[order,]
    fig=figure(figsize=(2,8))
    axmatrix=fig.add_axes([0.18,0.1,0.65,0.8])
    plt.register_cmap(name='cust_cmap', data=cdict)
    im=axmatrix.matshow(data,aspect='auto',origin='lower',cmap=plt.get_cmap('Reds'))
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    
    # Y axis
    aylabel=fig.add_axes([0.18,0.1,0.0,0.8])
    aylabel.set_yticks(yticks)
    aylabel.set_xticks([])
    
    # x axis
    axlabel=fig.add_axes([0.18,0.08,0.65,0.0])
    axlabel.set_xticks(xticks)
    axlabel.set_yticks([])
    
    #plot colorbar
    if showColBar:
        axcolor=fig.add_axes([0.86,0.1,0.03,0.8])
        colorbar(im,cax=axcolor)
    
    return fig

def ReadIntervals(bedfile):
    intervals =[]
    for tbed in IO.ColumnReader(sys.argv[1],ftype="bed"):
        intervals.append(tbed)
    return intervals
    
def ReadBigWig(intervals,wigfile):
    data = np.zeros((len(intervals),intervals[0].length()),dtype=float)
    for i,tbed in enumerate(intervals):
        data[i,:] += tbed.getWig(sys.argv[2],byPos=True)
    return data
    
# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" *.bed *.bw *.pdf/png")
    # Read intervals
    intervals = ReadIntervals(sys.argv[1])
    data = ReadBigWig(intervals,sys.argv[2])
    
    # process data
    data = standarize(data)
    # plot data
    fig=Plot(data,[-2000,-1000,0],[])
    
    # Save figure
    fig.savefig(sys.argv[3],format=sys.argv[3].split(".")[1].upper())
    


