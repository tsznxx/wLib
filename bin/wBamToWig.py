#!/usr/bin/python
#Last-modified: 06 Jan 2014 04:08:03 PM

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
import pysam
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
    p=argparse.ArgumentParser(description='Convert BAM file to Wiggle file in specific binsize. Reads are counted multiple times if not uniquely mapped. Contact Yunfei Wang to report any bugs (yfwang0405@gmail.com).',epilog='dependency pysam')
    p._optionals.title = "Options"
    p.add_argument("-t",dest="treatments",type=str,metavar="chip.bam", required=True, nargs='+',help="The ChIP experiment BAM file(s).")
    p.add_argument("-c",dest="control",type=str,metavar="control",default=None,help="The ChIP control BAM file.")
    p.add_argument("-b",dest="binsize",type=int,metavar='binsize',default=50,help="Binsize  used to generate the Wiggle file. [default=50]")
    p.add_argument("-e",dest="extend", type=int,metavar='extend' ,default=150,help="Extend read length by strand. Set \"0\" if no extension. [default= 150bp]")
    p.add_argument("-n",dest="normedto",type=int,metavar="normedto",default=10,help="Reads normalization . Set \"0\" if no normalization. Set \"10\" if want to normalize to 10 M reads. [default= 10]")
    p.add_argument("-R","--RNASeq",dest="RNASeq",action='store_true',default=False,help="Input is RNA Seq reads. No extension. [default= False].")
    p.add_argument("-P","--Paired",dest="Pairend",action='store_true',default=False,help="Input is paired end reads. No extension. [default= False].")
    p.add_argument("-s","--split",dest="split",action='store_true',default=False,help="Split wiggle into separate files by chromosomes. [default= False].")
    p.add_argument("-f","--forcestrand",dest="forcestrand",action='store_true',default=False,help="Force strand. [default= False].")
    p.add_argument("-p",dest="prefix",type=str,metavar="prefix",default="",help="Prefix to the output file name in format: Prefix_BamFileName_CtlBamFileName.wig. [default= \"\"]")
    p.add_argument("-o",dest="outputdir",type=str,metavar="ourputdir",default=".",help="Directory to put wiggle file(s). [default= \".\"].")
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

def bamToBins(bname,binsize=50,extend=150,normedto=10000000,paired=False,forcestrand=False):
    ''' Extend BAM and pileup to bins. '''
    # strand check
    if forcestrand:
        bins={"+":{},"-":{}}
    else:
        bins={".":{}}
    # Open bamfile
    sam=pysam.Samfile(bname,'rb')
    print >>sys.stderr, "BAM->bins ....."
    starttime=time.clock()
    ratio= normedto==0 and 1.0 or float(normedto)/sam.mapped
    # Bins initiation
    for tid,chrlen in enumerate(sam.lengths):
        for strand in bins:
            bins[strand][sam.references[tid]]=numpy.zeros(chrlen/binsize+extend/binsize+2,'float') # the number of bins > chromlen/binsize to simplify computation
    
    # For single end: Extension
    if not paired:
        cnt=0
        start=stop=0
        strand="."
        for read in sam: # sam is 0-based
            if not read.is_unmapped:
                cnt+=1
                if extend:
                    if read.is_reverse:
                        start, stop = max(0,read.pos+read.qlen-extend), read.pos+read.qlen-1
                    else:
                        start, stop = read.pos, read.pos+extend-1
                else:
                    start, stop = read.pos, read.pos+read.qlen-1
                # Add read to bins
                chrom, startbin, stopbin= sam.references[read.tid], start/binsize, stop/binsize
                strand = forcestrand and (read.is_reverse and "-" or "+") or "."
                # For startbin: 
                # take base pos 6 (1 based) as an example, binsize=50, readlen=60
                # pos 6-50 are in the startbin (50-6%50+1)
                # BAM file are 0 based, so bam.pos=5
                # So the bases in the startbin is binsize-bam.pos%binsize =
                bins[strand][chrom][startbin]-=start%binsize 
                # For stopbin:
                # stop is pos 65 (1 based, 6+60-1). bam.stop is 64 (0 based)
                # pos 51-65  are in the stopbin ((65-1)%50+1=15) Note: it is different from 65%50=15.
                # So the bases in the stopbin is bam.stop%binsize+1 = 64%50+1
                bins[strand][chrom][stopbin] +=stop%binsize+1 
                bins[strand][chrom][startbin:stopbin]+=binsize
        print >>sys.stderr, "Prcessed ",cnt, "reads.     "
    # For paired end: consider overlap
    else:
        cnt=0
        reads={}
        for read in sam:
            cnt+=1
            if read.is_unmapped==False and read.mate_is_unmapped==False and read.tid==read.rnext: # both mapped and in the same chrom
                if reads.has_key(read.qname): # has pair in reads
                    # find the mate pair
                    materead=None
                    for idx,tread in enumerate(reads[read.qname]):
                        if tread.tid == read.tid and tread.pos == read.pnext: # same chrom, mated
                            materead=tread
                            mateidx=idx
                    if materead: # if mated
                        # check the strand and pair
#                       if (read.is_read1 and materead.is_read2) or (read.is_read2 and materead.is_read1):
#                           print "Paired"
#                       else:
#                           print "not paired"
#                       if (read.is_reverse!=materead.is_reverse):
#                           print "complementary"
#                       else:
#                           print "same strand"

                        strand= forcestrand and ((read.is_read1 != read.is_reverse) and "+" or "-") or "."
                        #processing CIGAR  17M204733N83M=[(0, 17), (3, 204733), (0, 83)]
                        beds=[]
                        # one read
                        start=read.pos
                        for mtype,mlen in read.cigar:
                            if mtype==0: # match
                                beds.append((start,start+mlen-1))
                            elif mtype==1: # insertion
                                start-=mlen
                            start+=mlen
                        # the other read
                        start=materead.pos
                        for mtype,mlen in materead.cigar:
                            if mtype==0: # match
                                beds.append((start,start+mlen-1))
                            elif mtype==1: # insertion
                                start-=mlen
                            start+=mlen
                        # Add bed to bins
                        for start,stop in Merge(beds): # the beds are merged to avoid overlap between paired ends.
                            chrom, startbin, stopbin= sam.references[read.tid], start/binsize, stop/binsize
                            try:
                                bins[strand][chrom][startbin]-=start%binsize
                                bins[strand][chrom][stopbin] +=stop%binsize+1
                                bins[strand][chrom][startbin:stopbin]+=binsize
                            except:
                                print read
                                print reads[read.qname]
                                print read.cigar,",",reads[read.qname].cigar
                                print chrom, start, stop, startbin, stopbin
                                sys.exit()
                        # Delete read from reads
                        del reads[read.qname][mateidx]
                        if len(reads[read.qname])==0:
                            del reads[read.qname]
                    else: # not mate
                        reads[read.qname].append(read)
                else:
                    reads[read.qname]=[read]
        print >>sys.stderr, "Prcessed ",cnt, "reads.     "

    # Normalization
    print >>sys.stderr, "Normalization ratio=%-3.2f" % ratio
    for strand in bins:
        for chrom in bins[strand]:
            bins[strand][chrom]*=float(ratio)/binsize
            if sam.lengths[sam.references.index(chrom)]%binsize: # Set the last incomplete bin to 0
                bins[strand][chrom][-1]=0
    sam.close()

    # time used
    endtime=time.clock()
    print >>sys.stderr, "Time used: %-3.2f secs." % (endtime-starttime)
    return bins

def binsToWig(prefix,outputdir,samfile,bins,control,ctlbins=None,chrsplit=False):
    ''' Print to wig files. '''
    starttime=time.clock()
    # Get bname prefix
    bname=os.path.basename(samfile)
    bname=os.path.splitext(bname)[0]
    if prefix!="":
        bname=prefix+"_"+bname
    bname=outputdir+"/"+bname
    if control:
        bname += "_"+os.path.splitext(os.path.basename(control))[0]
    
    # chroms
    chroms={}
    sam=pysam.Samfile(samfile,'rb')
    for tid,length in enumerate(sam.lengths):
        chroms[sam.references[tid]]=length
    sam.close()

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
            bins[strand][chrom][length/binsize:] = 0 # ignore the bins >length/binsize
            for i in numpy.nonzero(bins[strand][chrom])[0]:
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
    
    # read control
    ctlbins=None
    if args.control:
        print >>sys.stderr, "Reading control:",args.control
        ctlbins=bamToBins(args.control,binsize,extend,normedto,isPaired,forcestrand)
    
    # read treatment
    for samfile in treatments:
        print >>sys.stderr, "Reading treatment:",samfile
        # Bam to Bins
        treatbins=bamToBins(samfile,binsize,extend,normedto,isPaired,forcestrand)
        # Wiggle output
        binsToWig(prefix,outputdir,samfile,treatbins,args.control,ctlbins,chromsplit)

        


