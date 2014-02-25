#!/usr/bin/python
#Last-modified: 29 Oct 2013 11:19:53 AM

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
import argparse
from wLib import Bed,IO,Utils

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Extend Bed/GeneBed region to upstream and/or downstream.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='ifname',type=str,metavar="input.bed",required=True,help="Input file. Can be stdin.")
    p.add_argument("-f","--format",dest="ftype",type=str,metavar="bed",default="bed",help="Format of input file. Default is 'bed'. Can be 'bed3', 'bedgraph', 'bed','peak','wig', 'sam2bed' or 'genepred'.")
    p.add_argument("-g","--genome",dest='genome',type=str,metavar="Genome",default=None,help="Genome version (hg19, mm10 .etc) or genome size file with chrom and size in each line.")
    p.add_argument("-u","--up",dest="up",type=int,metavar="upstream",default=0,help="bps extended to upstream. If minus, trim the 5' end.")
    p.add_argument("-d","--down",dest="down",type=int,metavar="downstream",default=0,help="bps extended to downstream. If minus, trim the 3' end.")
    p.add_argument("-o","--output",dest="ofname",type=str,metavar="output.bed",default="stdout",help="Output file. Default is stdout.")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    # Get parameters
    args=argParser()
    fh = IO.mopen(args.ofname, 'w')
    if args.genome:
        genome=Utils.genomeSize(args.genome)
    for item in IO.BioReader(args.ifname,args.ftype):
        tbed=item.extend(args.up,args.down, args.genome)
        print >> fh, tbed
    IO.mclose(fh)
