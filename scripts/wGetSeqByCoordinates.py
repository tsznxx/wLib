#!/usr/bin/python
#Last-modified: 29 Oct 2013 11:13:41 AM

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

import re
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
    p=argparse.ArgumentParser(description='Get a fragment from a Fasta file.',epilog='dependency wLib')
    p._optionals.title = "Options"
    p.add_argument("-i","--input",dest='ifname',type=str,metavar="input.fa",required=True,help="Fasta file name.")
    p.add_argument("-r","--region",dest='region',type=str,metavar="chr1:100-200:+",required=False,help = 'Chromosome region. Leave it empty if not applicable, i.e. "chr1:100-:-".')
    p.add_argument("-c","--chrom",dest='chrom',type=str,metavar="chrom",required=False,help="chromosome name.")
    p.add_argument("-s","--start",dest='start',type=int,metavar="start",required=False,default=None,help="start coordinate. Default: begining of the chromosome.")
    p.add_argument("-e","--end",dest='end',type=int,metavar="end",required=False,default=None,help="end coordiante. Default: end of the chromosome.")
    p.add_argument("-t","--strand",dest='strand',type=str,metavar='strand',required=False,default="+",help='strand. Default: "+"')
    p.add_argument("-o","--output",dest='ofname',type=str, default="stdout", metavar="output.fa",required=False,help="Output file name. Default: stdout")
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
    # Parse chromosome region
    if args.region:
        m = re.search("(\S+):(\S*)-(\S*):(\S*)",args.region)
        chrom, start, end, strand = m.groups()
    else:
        chrom = args.chrom
        start = args.start
        end   = args.end
        strand= args.strand
    if not chrom:
        raise ValueError("chromosome name is required.")
    start = int(start) if start else None
    end  = int(end) if end else None
    strand = strand if strand in ["+","-","."] else "+"
    # Get seq from Fasta file
    fa = FastaFile(args.ifname)
    seq = fa.getSeq(chrom,start,end,strand)
    fa.close()
    # Print result
    print >>fh, ">{0}:{1}-{2}:{3}".format(chrom, start and start or "", end and end or "", strand)
    print >>fh, seq.formatSeq()
    IO.mclose(fh)
