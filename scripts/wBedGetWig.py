#!/usr/bin/python
#Last-modified: 07 Mar 2013 01:37:19 PM

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

import sys
import string
from wLib.wBed import IO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        sys.exit("Example:"+sys.argv[0]+" *.bed *.bw ")
    for item in IO.ColumnReader(sys.argv[1],ftype='bed'):
        wigs=item.getWig(sys.argv[2])
        if len(wigs)>0:
            w=[wigs[i][2] for i in xrange(len(wigs))]
            item.score=max(w)
            bestpos=w.index(item.score)
        print item+"\t"+str(bestpos)

