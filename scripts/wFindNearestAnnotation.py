#!/usr/bin/python
#Last-modified: 29 Sep 2013 10:11:26 PM

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
from wLib import IO,GeneBedList,BedList,Utils

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
        print "Usage: "+sys.argv[0]+" annotation.tab/bed genomesize *.bed"
        print "       Find the nearest annotation for given bed."
    else:
        # check file
        if '.tab' in sys.argv[1]:
            ftype='gene'
        else:
            ftype='bed'

        # initiation annotations.
        annos={}
        #for chrom in IO.genomeSize('hg19'):
        for chrom in Utils.genomeSize(sys.argv[2]):
            if ftype=='bed':
                annos[chrom]=BedList()
            else:
                annos[chrom]=GeneBedList()
        
        # read annotations.
        for anno in IO.BioReader(sys.argv[1],ftype=ftype):
            if annos.has_key(anno.chrom):
                annos[anno.chrom].append(anno)

        # sort
        for chrom in annos:
            annos[chrom].sort()

        # Find nearest annoations
        for item in IO.BioReader(sys.argv[3],ftype='bed'):
            if annos.has_key(item.chrom) and len(annos[item.chrom])>0:
                tanno=annos[item.chrom][annos[item.chrom].bisect(item)]
                olen=item.overlapLength(tanno)
                if olen==1:
                    annostr=tanno.id+':overlap'
                    olen=0
                elif olen>1:
                    annostr=tanno.id+':body'
                elif item.stop>tanno.start:
                    if tanno.strand=="+":
                        annostr=tanno.id+':upstream'
                    else:
                        annostr=tanno.id+':downstream'
                        olen=-olen
                else:
                    if tanno.strand=="+":
                        annostr=tanno.id+':downstream'
                        olen=-olen
                    else:
                        annostr=tanno.id+':upstream'
                print str(item)+"\t"+str(tanno)+"\t"+str(olen)+"\t"+annostr
    
