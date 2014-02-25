#!/usr/bin/python
#Last-modified: 29 Oct 2013 11:00:37 AM

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
import argparse
from wLib import IO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Calculate GC content of Fasta file.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='fname',type=str,metavar="input.fa",required=True,help='Fasta file name. Can be "stdin".')
    p.add_argument("-o","--output",dest='ofname',type=str, default="stdout", metavar="output.gc",required=False,help="GC content file name. Default = stdout.")
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
    tlen=0
    gccnt=0
    # Read fasta file
    for tseq in IO.BioReader(args.fname,'fasta'):
        seq=str(tseq.seq).upper()
        cnt=seq.count("G")+seq.count("C")
        length=len(seq)-seq.count("N")
        print >> fh, "%s\t%-3.3f" % (tseq.id, float(cnt)/(length and length or 1))
        gccnt+=cnt
        tlen+=length
    print >> fh, "Total\t%-3.3f" % (float(gccnt)/(tlen and tlen or 1))
    IO.mclose(fh)
