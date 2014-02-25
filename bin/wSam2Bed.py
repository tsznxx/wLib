#!/usr/bin/python
#Last-modified: 20 Apr 2013 05:56:53 PM

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
import pysam

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
        sys.exit("Example:"+sys.argv[0]+" *.bam/sam >*.bed ")
    sam=pysam.Samfile(sys.argv[1], sys.argv[1].endswith("bam") and 'rb' or 'r')
    for read in sam:
        if not read.is_unmapped:
            print "%s\t%d\t%d\t%s\t%d\t%s" %(sam.references[read.tid],read.pos,read.aend,read.qname,read.mapq,read.is_reverse and "-" or "+")
    sam.close()

