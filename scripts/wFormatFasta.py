#!/usr/bin/python
#Last-modified: 29 Oct 2013 10:56:41 AM

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
from wLib import Utils,IO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Format sequences in Fasta file to fixed length.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='ifname',type=str,metavar="input.fa",required=True,help='Fasta file name. Can be "stdin".')
    p.add_argument("-l","--length",dest='length',type=int,metavar="length",default=100,help="Length of each line. Default is 100.")
    p.add_argument("-o","--output",dest='ofname',type=str, default="stdout", metavar="output.fa",required=False,help="Output file name. Default is stdout.")
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
    args = argParser()
    fh = IO.mopen(args.ofname,'w')
    for fa in IO.BioReader(args.ifname,'fasta'):
        print >>fh, ">"+fa.id
        print >>fh, fa.seq.formatSeq(args.length)
    IO.mclose(fh)
