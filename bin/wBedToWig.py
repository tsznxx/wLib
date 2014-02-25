#!/usr/bin/python
#Last-modified: 05 Nov 2013 02:24:52 PM

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

import os
import sys
import string
import argparse
import numpy
import time

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Convert Bed file to Wiggle file in specific binsize. This program will multiply the score of each bed. There are bugs if allowed mutilple hits in mapping.')
    p._optionals.title = "Options"
    p.add_argument("-t",dest="treatments",type=str,metavar="sample1.bed", required=True, nargs='+',help="The experiment Bed file(s).")
    p.add_argument("-c",dest="control",type=str,metavar="control",help="The control Bed file.")
    p.add_argument("-g",type=str,metavar='genome',required=True,help="Chrom sizes file.")
    p.add_argument("-b",dest="binsize",type=int,metavar='binsize',default=50,help="Binsize  used to generate the Wiggle file. [default=50]")
    p.add_argument("-e",dest="extend", type=int,metavar='extend' ,default=150,help="Extend read length by strand. Set \"0\" if no extension. [default= 150bp]")
    p.add_argument("-n",dest="normedto",type=int,metavar="normedto",default=10,help="Reads normalization . Set \"0\" if no normalization. Set \"10\" if want to normalized to 10 M reads. [default= 10 M]")
    p.add_argument("-R","--RNASeq",dest="RNASeq",action='store_true',default=False,help="Input is RNA Seq reads. No extension by default.")
    p.add_argument("-P","--Paired",dest="Pairend",action='store_true',default=False,help="Input is paired end reads. No extension by default.")
    p.add_argument("-s","--split",dest="split",action='store_true',default=False,help="Split wiggle into separate files by chromosomes.")
    p.add_argument("-f","--forcestrand",dest="forcestrand",action='store_true',default=False,help="Force strand.")
    p.add_argument("-p",dest="prefix",type=str,metavar="prefix",default="",help="Prefix to the output file name in format: Prefix_BedFileName.wig. Default is BedFileName.wig.")
    p.add_argument("-o",dest="outputdir",type=str,metavar="ourputdir",default=".",help="Directory to put wiggle file(s). Default is current directory.")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def Merge(beds):
    newbeds=[]
    newbeds.append(list(beds[0]))
    for start,stop in beds[1:]:
        if newbeds[-1][0]<= start <=newbeds[-1][1]:
            newbeds[-1][1]=max(newbeds[-1][1],stop)
        else:
            newbeds.append([start,stop])
    return newbeds

def bamToBins(bname,chroms,binsize=50,extend=150,normedto=10000000,paired=False,forcestrand=False):
    ''' Extend BAM and pileup to bins. '''
    # strand check
    if forcestrand:
        bins={"+":{},"-":{}}
    else:
        bins={".":{}}
    # Open bamfile
    # sam=pysam.Samfile(bname,'rb')
    if bname == "stdin" or bname == "-":
        sam = sys.stdin
    else:
        sam = open(bname)
    print >>sys.stderr, "BAM->bins ....."
    starttime=time.clock()
    ratio= 1.0 #normedto==0 and 1.0 or float(normedto)/sam.mapped
    # Bins initiation
    for chrom in chroms:
        for strand in bins:
            bins[strand][chrom]=numpy.zeros(chroms[chrom]/binsize+extend/binsize+1,'float')
    # For single end: Extension
    if not paired:
        cnt=0
        start=stop=0
        strand="."
        for read in sam: # sam is 0-based
            cnt += 1
            read = read.split("\t")
            try:
                rscore = float(read[4])
            except:
                rscore =1.0
            if rscore <= 0: # we don't count beds with score <= 0
                continue
            chrom = read[0]
            rstart = int(read[1])
            rstop = int (read[2])
            try:
                strand = read[5]
            except:
                strand = "."
            if cnt%100000==0:
                print >>sys.stderr, "Processed ",cnt,"reads.      \r",
            if extend:
                if strand == "-":
                    start, stop = max(0,rstop-extend), rstop-1
                else:
                    start, stop = rstart, rstart+extend-1
            else:
                start, stop = rstart, rstop-1
            # Add read to bins
            startbin, stopbin= start/binsize, stop/binsize
            # For startbin: 
            # take base pos 6 (1 based) as an example, binsize=50, readlen=60
            # pos 6-50 are in the startbin (50-6%50+1)
            # BAM file are 0 based, so bam.pos=5
            # So the bases in the startbin is binsize-bam.pos%binsize =
            bins[strand][chrom][startbin]-=start%binsize*rscore
            # For stopbin:
            # stop is pos 65 (1 based, 6+60-1). bam.stop is 64 (0 based)
            # pos 51-65  are in the stopbin ((65-1)%50+1=15) Note: it is different from 65%50=15.
            # So the bases in the stopbin is bam.stop%binsize+1 = 64%50+1
            bins[strand][chrom][stopbin] +=(stop%binsize+1)*rscore
            bins[strand][chrom][startbin:stopbin]+=binsize*rscore
        print >>sys.stderr, "Prcessed ",cnt, "reads.     "
    # For paired end: consider overlap
#    else:
#        cnt=0
#        reads={}
#        for read in sam:
#            cnt+=1
#            if cnt%10000==0:
#                print >>sys.stderr, "Processed ",cnt,"reads.      \r",
#            if read.is_unmapped==False and read.mate_is_unmapped==False and read.tid==read.rnext: # both mapped and in the same chrom
#                if reads.has_key(read.qname): # has pair in reads
#                    # find the mate pair
#                    materead=None
#                    for idx,tread in enumerate(reads[read.qname]):
#                        if tread.tid == read.tid and tread.pos == read.pnext: # same chrom, mated
#                            materead=tread
#                            mateidx=idx
#                    if materead: # if mated
#                        # check the strand and pair
##                       if (read.is_read1 and materead.is_read2) or (read.is_read2 and materead.is_read1):
##                           print "Paired"
##                       else:
##                           print "not paired"
##                       if (read.is_reverse!=materead.is_reverse):
##                           print "complementary"
##                       else:
##                           print "same strand"
#
#                        strand= forcestrand and ((read.is_read1 != read.is_reverse) and "+" or "-") or "."
#                        #processing CIGAR  17M204733N83M=[(0, 17), (3, 204733), (0, 83)]
#                        beds=[]
#                        # one read
#                        start=read.pos
#                        for mtype,mlen in read.cigar:
#                            if mtype==0: # match
#                                beds.append((start,start+mlen-1))
#                            elif mtype==1: # insertion
#                                start-=mlen
#                            start+=mlen
#                        # the other read
#                        start=materead.pos
#                        for mtype,mlen in materead.cigar:
#                            if mtype==0: # match
#                                beds.append((start,start+mlen-1))
#                            elif mtype==1: # insertion
#                                start-=mlen
#                            start+=mlen
#                        # Add bed to bins
#                        for start,stop in Merge(beds): # the beds are merged to avoid overlap between paired ends.
#                            chrom, startbin, stopbin= sam.references[read.tid], start/binsize, stop/binsize
#                            try:
#                                bins[strand][chrom][startbin]-=start%binsize
#                                bins[strand][chrom][stopbin] +=stop%binsize+1
#                                bins[strand][chrom][startbin:stopbin]+=binsize
#                            except:
#                                print read
#                                print reads[read.qname]
#                                print read.cigar,",",reads[read.qname].cigar
#                                print chrom, start, stop, startbin, stopbin
#                                sys.exit()
#                        # Delete read from reads
#                        del reads[read.qname][mateidx]
#                        if len(reads[read.qname])==0:
#                            del reads[read.qname]
#                    else: # not mate
#                        reads[read.qname].append(read)
#                else:
#                    reads[read.qname]=[read]
#        print >>sys.stderr, "Prcessed ",cnt, "reads.     "

    # Normalization
    ratio= normedto==0 and 1.0 or float(normedto)/cnt
    print >>sys.stderr, "Normalization ratio=%-3.2f" % ratio
    for strand in bins:
        for chrom in bins[strand]:
            bins[strand][chrom]*=float(ratio)/binsize
            if chroms[chrom]%binsize: # Set the last incomplete bin to 0
                bins[strand][chrom][-1]=0
    if sam != sys.stdin:
        sam.close()

    # time used
    endtime=time.clock()
    print >>sys.stderr, "Time used: %-3.2f secs." % (endtime-starttime)
    return bins

def binsToWig(prefix,outputdir,chroms,bins,ctlbins=None,chrsplit=False):
    ''' Print to wig files. '''
    starttime=time.clock()
    # Get bname prefix
    bname=os.path.basename(samfile)
    bname=os.path.splitext(bname)[0]
    if prefix!="":
        bname=prefix+"_"+bname
    bname=outputdir+"/"+bname
    
    # binsize
    chrom=chroms.keys()[0]
    strand=bins.keys()[0]
    binsize=chroms[chrom]/len(bins[strand][chrom])+1
    
    # strand
    for strand in bins:
        if strand!=".":
            suffix= strand=="+" and "_plus.wig" or "_minus.wig"
        else:
            suffix=".wig"
        # Open wiggle file
        if not chrsplit:
            fh=open(bname+suffix,'w')
            print >>fh, "track type=wiggle_0 name="+os.path.basename(bname)
        
        # Print wiggle
        for chrom in sorted(chroms.keys()):
            length=chroms[chrom]
            print >>sys.stderr, "Bins->Wiggle ", chrom, strand, '  \r',
            
            # Open Wiggle file
            if chrsplit:
                fh=open(bname+"_"+chrom+suffix,'w')
                print >>fh, "track type=wiggle_0 name="+os.path.basename(bname)+"_"+chrom
            
            # normalized to control
            if ctlbins and ctlbins[strand].has_key(chrom):
                bins[strand][chrom]-=ctlbins[strand][chrom]
                bins[strand][chrom]+=abs(bins[strand][chrom])
                bins[strand][chrom]/=2
            
            # Print wiggle lines
            print >>fh, "variableStep chrom="+chrom+" span="+str(binsize)
            for i in numpy.nonzero(bins[strand][chrom])[0]: # xrange(length/binsize):
                print >>fh, "%d\t%-3.2f" % (i*binsize+1, bins[strand][chrom][i])
            
            # Close file handle
            if chrsplit:
                fh.close()
    
    # Close file handle
    if not chrsplit:
        fh.close()

    print >>sys.stderr
    endtime=time.clock()
    print >>sys.stderr, "Time used: %-3.2f secs." % (endtime-starttime)
    return

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    # Get parameters
    args=argParser()
    binsize=args.binsize
    genome=args.genome
    extend=args.extend
    normedto=args.normedto*1000000
    chromsplit=args.split
    treatments=args.treatments
    outputdir=args.outputdir
    prefix=args.prefix
    forcestrand=args.forcestrand
    isPaired=args.Pairend
    isRNASeq = args.RNASeq
    if isPaired or isRNASeq: # No extension if RNA Seq
        extend=0
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    
    # Print parameters
    print >>sys.stderr, "Parameters:"
    print >>sys.stderr, "\tbinsize=", binsize,"bp"
    print >>sys.stderr, isPaired and "\tPaired end reads." or "\tSingle end reads."
    print >>sys.stderr, extend   and "\textend to %d bp" % extend or "\tNo extension"
    print >>sys.stderr, normedto and "\tnormalized to %dM reads" % args.normedto or "\tNo normalization"
    if isRNASeq:
        print >>sys.stderr, "\tRNA Seq reads."
    print >>sys.stderr, forcestrand and "\tForce strand." or "\tDon't force strand."
    if chromsplit:
        print >>sys.stderr, "\tSplit wiggle by chromosomes."
    print >>sys.stderr, outputdir=="." and "\toutput in current dir" or "\toutput in "+outputdir
    print >>sys.stderr

    # read genome
    chroms={}
    with open(args.genome) as fh:
        for line in fh:
            chrom,size = line.split("\t")
            chroms[chrom]=int(size)

    # read control
    ctlbins=None
    if args.control:
        print >>sys.stderr, "Reading control:",args.control
        ctlbins=bamToBins(args.control,chroms,binsize,extend,normedto,isPaired,forcestrand)
    
    # read treatment
    for samfile in treatments:
        print >>sys.stderr, "Reading treatment:",samfile
        # Bam to Bins
        treatbins=bamToBins(samfile,chroms,binsize,extend,normedto,isPaired,forcestrand)
        # Wiggle output
        binsToWig(prefix,outputdir,chroms,treatbins,ctlbins,chromsplit)

        


