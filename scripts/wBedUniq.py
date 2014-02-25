#!/usr/bin/python
#Last-modified: 29 Oct 2013 01:45:29 PM

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

import sys,copy
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
    p=argparse.ArgumentParser(description='Get unique Beds from a Bed file. Bed file should be sorted first.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='ifname',type=str,metavar="sorted.bed",required=True,help="Input file. Should be sorted. Can be 'stdin'.")
    p.add_argument("-f","--format",dest="ftype",type=str,metavar="bed",default="bed",help="Format of input file. Default is 'bed'. Can be 'bed3', 'bedgraph', 'bed','peak','wig', 'sam2bed' or 'genepred'.")
    p.add_argument("-n","--names",dest="names",choices=['p','c','n'],default='n',help="'p': prefix_num,'c': collapse, id1;id2..., 'n': number.")
    p.add_argument("-p","--prefix",dest="prefix",type=str,metavar="prefix",default="item",help="Prefix of Bed names. Requires -n parameter.")
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
    lastbed = Bed("chr1\t0\t0")
    cnt = 0 # count of unique items
    for item in IO.BioReader(args.ifname,args.ftype):
        if item.chrom == lastbed.chrom and item.start == lastbed.start and item.stop == lastbed.stop:
            if args.names == 'c':
                lastbed.id += ";"+item.id
        else:
            if lastbed.stop !=0:
                cnt += 1
                if args.names == 'n':
                    lastbed.id = str(cnt)
                elif args.names == 'p':
                    lastbed.id = args.prefix+"_"+str(cnt)
                print >> fh, lastbed

            lastbed = copy.deepcopy(item)
    if lastbed.stop != 0:
        cnt += 1
        if args.names == 'n':
            lastbed.id = str(cnt)
        elif args.names == 'p':
            lastbed.id = args.prefix+"_"+str(cnt)
        print >> fh, lastbed
    IO.mclose(fh)
