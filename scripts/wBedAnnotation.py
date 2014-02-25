#!/usr/bin/python
#Last-modified: 06 Feb 2013 03:41:23 PM

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
from wLib.wBed import IO,BedMap

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
        print "Usage: "+sys.argv[0]+" *.tab *.bed"
    else:
        # reading GeneTab file
        genes=BedMap('hg19')
        if 'tab' in sys.argv[1][-3:]:
            genes.loadBedToMap(sys.argv[1],'gene')
        else:
            genes.loadBedToMap(sys.argv[1],'bed')
        
        # reading bed file
        for bed in IO.ColumnReader(sys.argv[2],ftype='bed'):
            overlaps=genes.intersectBed(bed)
            anno=None
            overmax=0
            if overlaps:
                print str(bed)+"\t"+str(overlaps)
            else:
                print str(bed)+"\tNone\tNone\t0\t0\t0\t0"

                

