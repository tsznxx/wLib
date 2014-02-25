#!/usr/bin/python
#Last-modified: 29 Oct 2013 11:16:36 AM

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
from wLib import IO,FastaFile

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Get sequences by a list of names.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='ifname',type=str,metavar="input.fa",required=True,help='Fasta file name.')
    p.add_argument("-n","--names",dest='names',type=str,metavar="names.lst",required=True,help='A file with sequence names. Can be "stdin".')
    p.add_argument("-o","--output",dest='ofname',type=str, default="stdout", metavar="output.fa",required=False,help='Output file name. Default is "stdout".')
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
    fh = IO.mopen(args.ofname, 'w')
    # Fasta file
    fa = FastaFile(args.ifname)
    f = IO.mopen(args.names)
    for faid in f:
        faid=faid.rstrip()
        seq=fa.getSeq(faid)
        print >> fh, ">"+faid
        print >> fh, seq.formatSeq()
    fa.close()
    IO.mclose(f)
    IO.mclose(fh)

