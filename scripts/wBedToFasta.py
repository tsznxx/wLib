#!/usr/bin/python
#Last-modified: 29 Oct 2013 02:43:58 PM

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

import os,sys
import string
import argparse
from wLib import IO,Utils,FastaFile

# ------------------------------------
# constants
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Fetch sequences for Bed file from genome file.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='ifname',type=str,metavar="input.bed",required=True,help="Input file. Can be 'stdin'.")
    #p.add_argument("-s","--strand",dest="strand",type=str, metavar="strand",default=".",choices=".+-=",help=".: take the strand as indicated in input file. +: take the plus strand of the genome. -: take the minus sequence of the genome. =: take both strands.")
    p.add_argument("-g","--genome",dest="genome",type=str,metavar="genome",required=True,help="Genome file.")
    p.add_argument("-f","--format",dest="ftype",type=str,metavar="bed",default="bed",help="Format of input file. Default is 'bed'. Can be 'bed3', 'bedgraph', 'bed','peak','wig', 'sam2bed', 'any' or 'genepred'.")
    p.add_argument("-c","--case",dest="case",type=str,choices=['o', 'u','l'],default="o",help="o: original(default), u: uppercase, l: lowercase")
    p.add_argument("-l","--linelen",dest='linelength',type=int,metavar='linelength',default=100,help="Bases in each line. Default is 100. 0 for unlimited length.")
    p.add_argument("-o","--output",dest="ofname",type=str,metavar='output.fa',default="stdout",help="Output file. Default is 'stdout'.")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

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
    # Get parameters
    args=argParser()
    fh = IO.mopen(args.ofname)
    genome=FastaFile(args.genome)
    for i,item in enumerate(IO.BioReader(args.ifname,args.ftype)):
        try:
            strand = item.strand
        except:
            strand = "+"
        seq=genome.getSeq(item.chrom,item.start,item.stop,strand)
        if len(seq)>0:
            print >> fh, '>'+(item.id!="NONAME" and item.id or "item_"+str(i))
            if args.linelength:
                seq = seq.formatSeq(args.linelength)
            if args.case=='u':
                seq=seq.upper()
            elif args.case=='l':
                seq=seq.lower()
            print >> fh, seq
    genome.close()
    IO.mclose(fh)
