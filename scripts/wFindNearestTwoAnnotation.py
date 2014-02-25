#!/usr/bin/python
#Last-modified: 06 Feb 2013 03:43:03 PM

""" Module/Scripts Description

Copyright (c) 2008 Yunfei Wang <tszn1984@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: 1.0.0
@author:  Yunfei Wang
@contact: tszn1984@gmail.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import sys
import string
from wLib.wBed import IO
import bisect

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
        print "Usage: "+sys.argv[0]+" annotation.tab/bed *.bed"
        print "       Find the nearest annotation for given bed."
    else:
        # check file
        if sys.argv[1][-3:]=='bed':
            ftype='bed'
        elif sys.argv[1][-3:]=='tab':
            ftype='gene'
        else:
            print >>sys.stderr, "File type error: file type should be either '.bed' or '.tab'."
            sys.exit(0)

        # initiation annotations.
        annos={}
        for chrom in {"K12":4639675}:
             annos[chrom]=[]
        
        # read annotations.
        for anno in IO.ColumnReader(sys.argv[1],ftype):
            annos[anno.chr].append(anno)

        # sort
        for chrom in annos:
            annos[chrom].sort()

        # Find nearest annoations
        for item in IO.ColumnReader(sys.argv[2],'bed'):
            tanno=bisect.bisect(annos[item.chr],item)
            print str(item)+"\t"+annos[item.chr][tanno-1].id+"\t"+annos[item.chr][tanno-1].strand+"\t"+str(item.overlapLength(annos[item.chr][tanno-1]))+"\t"+annos[item.chr][tanno].id+"\t"+annos[item.chr][tanno].strand+"\t"+str(item.overlapLength(annos[item.chr][tanno]))

