#!/usr/Bin/python
#Last-modified: 14 Dec 2012 02:50:19 PM

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
import ghmm
from wLib import lprob
import wWigIO
from wLib.wBed import Utils,IO,Peak
from wLib.fdr import storey_qvalues


# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def argParser():
    ''' Parse arguments. '''
    p=argparse.ArgumentParser(description='Convert BAM file to Wiggle file in specific Binsize.',epilog='dependency pysam')
    p._optionals.title = "Options"
    p.add_argument("-t",dest="treatments",type=str,metavar="chip.bam", required=True, nargs='*',help="The ChIP experiment BAM file(s). Histone marks supported are \"H3K4me1\", \"H3K4me3\", \"H3K27ac\", \"H3K27me3\", \"H3K36me3\". Histone mark names should be shown in file names.")
    p.add_argument("-c",dest="control",type=str,metavar="control",help="The ChIP control BAM file.")
    p.add_argument("-b",dest="Binsize",type=int,metavar='Binsize',default=50,help="Binsize  used to generate the Wiggle file. [default=50]")
    p.add_argument("-e",dest="extend", type=int,metavar='extend' ,default=150,help="Extend read length by strand. Set 0 if no extension. [default= 150bp]")
    p.add_argument("-n",dest="normedto",type=int,metavar="normedto",default=1000000,help="Reads normalization . Set 0 if no normalization. [default= 1M]")
    p.add_argument("-s",dest="split",action='store_true',default=False,help="\"-w\" required. Split wiggle into separate files by chromosomes.")
    p.add_argument("-o",dest="outputdir",type=str,metavar="ourputdir",default=".",help="Directory to put wiggle file(s). Default is current directory.")
    p.add_argument("-w",dest="wig",action='store_true',default=False,help="Save wiggle file. May use together with \"-s\".")
    if len(sys.argv)==1:
        sys.exit(p.print_help())
    args = p.parse_args()
    return args

def printTime():
    ''' Print time in fixed format. '''
    return time.strftime("#### %d %b %Y %H:%M:%S ####", time.gmtime())

class BamBins:
    ''' Bam file to Bins. '''
    def __init__(self):
        ''' Initiate BamBins. '''
        self.Binsize=50

    def bamToBins(samfile,Binsize=50,extend=150,normedto=1000000):
        ''' Extend BAM and pileup to Bins. '''
        Bins={}
        sam=pysam.Samfile(samfile,'rb')
        print >>sys.stderr, printTime(), "BAM->Bins ....."
        ratio= normedto==0 and 1.0 or float(normedto)/sam.mapped
        for tid,chrlen in enumerate(sam.lengths):
            Bins[sam.references[tid]]=numpy.zeros(chrlen/Binsize+extend/Binsize+1,'float')
        
        # Extension
        cnt=0
        start=stop=0
        for read in sam: # sam is 0-based
            if not read.is_unmapped:
                cnt+=1
                if cnt%10000==0:
                    print >>sys.stderr, "Processed ",cnt,"reads.      \r",
                if extend:
                    if read.is_reverse:
                        start, stop = max(0,read.pos+read.qlen-extend), read.pos+read.qlen-1
                    else:
                        start, stop = read.pos, read.pos+extend-1
                else:
                    start, stop = read.pos, read.pos+read.qlen-1
                # Add read to Bins
                chrom, startBin, stopBin = sam.references[read.tid], start/Binsize, stop/Binsize
                # For startBin: 
                # take base pos 6 (1 based) as an example, Binsize=50, readlen=60
                # pos 6-50 are in the startBin (50-6%50+1)
                # BAM file are 0 based, so bam.pos=5
                # So the bases in the startBin is Binsize-bam.pos%Binsize =
                Bins[chrom][startBin]-=start%Binsize 
                # For stopBin:
                # stop is pos 65 (1 based, 6+60-1). bam.stop is 64 (0 based)
                # pos 51-65  are in the stopBin ((65-1)%50+1=15) Note: it is different from 65%50=15.
                # So the bases in the stopBin is bam.stop%Binsize+1 = 64%50+1
                Bins[chrom][stopBin] +=stop%Binsize+1 
                Bins[chrom][startBin:stopBin]+=Binsize
        print >>sys.stderr, printTime(), "Prcessed ",cnt, "reads.     "
        
        
        # Normalization
        print >>sys.stderr, printTime(), "Normalization to", normedto, "reads. ratio=%-3.2f" % ratio
        for chrom in Bins:
            Bins[chrom]*=float(ratio)/Binsize
            if sam.lengths[sam.references.index(chrom)]%Binsize: # Set the last incomplete Bin to 0
                Bins[chrom][-1]=0
        sam.close()
    
        # time used
        print >>sys.stderr, printTime(), "Bam->Bins finished." 
        print >>sys.stderr
        return Bins
    bamToBins=staticmethod(bamToBins)
    def normalizeBins(treatBins,ctlBins=None):
        ''' Normalize Bins against control. '''
        if ctlBins:
            print >>sys.stderr, printTime(), "Normalized to control Bins."
            for chrom in treatBins:
                treatBins[chrom]-=ctlBins[chrom]
                treatBins[chrom]+=abs(treatBins[chrom])
                treatBins[chrom]/=2
        print >>sys.stderr
        return
    normalizeBins=staticmethod(normalizeBins)
    def BinsToWig(outputdir,samfile,Bins,chrsplit=False):
        ''' Print to wig files. '''
        # Get bname prefix
        bname=os.path.basename(samfile)
        bname=outputdir+"/"+os.path.splitext(bname)[0]
        
        # Binsize
        Binsize=self.Binsize
        
        # Open wiggle file
        if not chrsplit:
            fh=open(bname+".wig",'w')
            print >>fh, "track type=wiggle_0 name="+os.path.basename(bname)
        
        # Print wiggle
        for chrom in Bins:
            print >>sys.stderr, printTime(), "Bins->Wiggle ", chrom
            
            # Open Wiggle file
            if chrsplit:
                fh=open(bname+"_"+chrom+".wig",'w')
                print >>fh, "track type=wiggle_0 name="+os.path.basename(bname)+"_"+chrom
            
            # Print wiggle lines
            print >>fh, "variableStep chrom="+chrom+" span="+str(Binsize)
            for i,d in enumerate(Bins[chrom]):
                depth=float("%.2f" % d)
                if depth>0:
                    print >>fh, "%d\t%-3.2f" % (i*Binsize+1,depth)
            
            # Close file handle
            if chrsplit:
                fh.close()
        
        # Close file handle
        if not chrsplit:
            fh.close()
    
        print >>sys.stderr, printTime(), "Bins->Wiggle finished."
        print >>sys.stderr
        return
    BinsToWig=staticmethod(BinsToWig)
    def BinsFromBigWig(bwfile,binsize):
        ''' Read bins from BigWig file.'''
        Bins={}
        chroms=Utils.getBigWigChroms(bwfile)
        for chrom in chroms:
            Bins[chrom]=numpy.zeros(chroms[chrom]/binsize)
            for start,stop,val in IO.BigWigReader(bwfile,chrom):
                Bins[chrom][start/binsize]=val
        return Bins
    BinsFromBigWig=staticmethod(BinsFromBigWig)
    def callPeakFromBins(Bins,pvalue=1e-5):
        ''' Call peak with poisson distribution. Return the hmmStates, lamda and threshold for each chromosome.'''
        print >>sys.stderr, printTime(), "Call peak for each chromosome."
        hmmStates={}
        lams={}
        for chrom in Bins:
            # All bin values are multiplied by 100, because the threshold is usually integers.
            lam = sum(Bins[chrom])/len(Bins[chrom])*100
            # find the optimal threshold by binary search
            lthre=0.0
            hthre=max(Bins[chrom])*100.0  # 
            while hthre-lthre>0.01:
                pval = lprob.poisson_cdf((hthre+lthre)/2,lam,False)
                if pval > pvalue:
                    lthre=(hthre+lthre)/2
                else:
                    hthre=(hthre+lthre)/2
            # Print threshold
            print >>sys.stderr, printTime(), chrom, "lambda=",lam,"threshold=",hthre
            # convert to binary value
            hmmStates[chrom]=(numpy.greater_equal(Bins[chrom],hthre/100)-0).tolist()
            lams[chrom]=lam
        print >>sys.stderr, printTime(), "Call peak finished."
        print >>sys.stderr
        return (hmmStates,lams)
    callPeakFromBins=staticmethod(callPeakFromBins)
    def trainHMM(hmmState):
        ''' Train HMM with the given chromosome. '''
        print >>sys.stderr, printTime(), "Train HMM with one chromosome."
        T=[[0.9,0.1],[0.1,0.9]]
        e1=[0.1,0.9]
        e0=[0.9,0.1]
        E=[e0,e1]
        pi=[0.9,0.1] # initial 10% are peak?
        sigma=ghmm.IntegerRange(0,2) # 0, 1
        m = ghmm.HMMFromMatrices(sigma,ghmm.DiscreteDistribution(sigma),T,E,pi)
        m.baumWelch(ghmm.EmissionSequence(sigma,hmmState))
        print >>sys.stderr, printTime(), "Train HMM finished."
        print >>sys.stderr
        return m
    trainHMM=staticmethod(trainHMM)
    def decodeHMM(m,hmmStates):
        ''' Decode HMM. '''
        print >>sys.stderr, printTime(), "Decode HMM for each chromosome."
        sigma=ghmm.IntegerRange(0,2)
        for chrom in hmmStates:
            print >>sys.stderr, printTime(), chrom
            state, score = m.viterbi(ghmm.EmissionSequence(sigma,hmmStates[chrom]))
            hmmStates[chrom]=state
        print >>sys.stderr, printTime(), "Decode HMM finished."
        print >>sys.stderr
    decodeHMM=staticmethod(decodeHMM)
    def parseHMM(hmmStates,Binsize):
        ''' Parse HMM to peaks. '''
        for chrom in hmmStates:
            print >>sys.stderr, printTime(), chrom
            flag=False
            for i,s in enumerate(hmmStates[chrom]):
                if s==1:  # current =1
                    if flag: # last =1 
                        stop+=1
                    else:    # last =0 ,start of peak
                        start=i
                        stop=i
                        flag=True
                else:    # curent =0
                    if flag:  # last =1, end of peak
                        yield (chrom,start*Binsize,(stop+1)*Binsize)
                        flag=False
            if flag: # end of chromosome
                yield (chrom,start*Binsize,(stop+1)*Binsize)
        return
    parseHMM=staticmethod(parseHMM)
    def HMMToPeaks(hmmStates,treatbins,Binsize):
        ''' Parse HMM result to Peaks. '''
        # parse HMM and generate peak files.
        print >>sys.stderr, printTime(), "Parse Peaks from HMM result."

        peaks=[]
        pval=[]
        cnt=1
        for chrom,start,stop in BamBins.parseHMM(hmmStates,Binsize):
            tpeak=Peak([chrom,start,stop,histmark+"_"+str(cnt)])
            bins=treatBins[chrom][start/Binsize:stop/Binsize]
            tpeak.score=max(bins)
            tpeak.signalvalue=numpy.mean(bins)
            tpeak.pvalue = -10*numpy.log10(lprob.poisson_cdf(tpeak.score*100,lams[chrom],False))
            tpeak.peak=bins.argmax()*Binsize+Binsize/2
            # append to peaks
            peaks.append(tpeak)
            #pval.append(tpeak.pvalue)
            cnt+=1
        #qval=storey_qvalues(pval)
        #for i in xrange(len(peaks)):
            #peaks[i].pvalue=-10*numpy.log10(pval[i])
            #peaks[i].qvalue=-10*numpy.log10(qval[i])
        print >>sys.stderr, printTime(), "Parse Peaks finished."
        print >>sys.stderr
        return peaks
    HMMToPeaks=staticmethod(HMMToPeaks)
    

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    # Get parameters
    args=argParser()
    Binsize=args.Binsize
    extend=args.extend
    normedto=args.normedto
    chromsplit=args.split
    outputdir=args.outputdir
    towig=args.wig
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)
    histmarks="H3K4me3,H3K4me1,H3K27ac,H3K27me3,H3K36me3".split(",")
    treatments={}
    for f in args.treatments:
        for h in histmarks:
            if h.upper() in f.upper():
                treatments[h]=f
    
    # Print parameters
    print >>sys.stderr, printTime(), "Parameters:"
    print >>sys.stderr, printTime(), "    Binsize=", Binsize,"bp"
    print >>sys.stderr, printTime(), extend   and "    extend to %d bp" % extend or "    No extension"
    print >>sys.stderr, printTime(), normedto and "    normalized to %d reads" % normedto or "    No normalization"
    if chromsplit:
        print >>sys.stderr, printTime(), "    Split wiggle by chromosomes."
    print >>sys.stderr, printTime(), outputdir=="." and "    output in current dir" or "    output in "+outputdir
    print >>sys.stderr
    
    # read control
    ctlBins=None
    if args.control:
        print >>sys.stderr, printTime(), "Reading control:",args.control
        if args.control.endswith("bam"):
            print >>sys.stderr, printTime(), "File format is BAM."
            ctlBins=BamBins.bamToBins(args.control,Binsize,extend,normedto)
        elif args.control.endswith("bw") or args.control.endswith("bigwig"):
            print >>sys.stderr, printTime(), "File format is BigWig."
            ctlBins=BamBins.BinsFromBigWig(args.control,Binsize)
        else:
            sys.exit(printTime()+ "ERROR: "+args.control+" file format not suported. Options are *.bam, *.bw or *.bigwig.")            
    
    # read treatment
    Peaks={}
    for histmark in treatments:
        print >>sys.stderr, printTime(), "Reading treatment:",treatments[histmark]
        # BAM to Bins
        if treatments[histmark].endswith("bam"):
            print >>sys.stderr, printTime(), "File format is BAM."
            treatBins=BamBins.bamToBins(treatments[histmark],Binsize,extend,normedto)
        elif treatments[histmark].endswith("bw") or treatments[histmark].endswith("bigwig"):
            print >>sys.stderr, printTime(), "File format is BigWig."
            treatBins=BamBins.BinsFromBigWig(treatments[histmark],Binsize)
        else:
            sys.exit(printTime()+ "ERROR"+treatments[histmark]+" file format not suported. Options are *.bam, *.bw or *.bigwig.")
        
        # normalize Bins
        BamBins.normalizeBins(treatBins,ctlBins)
        
        # Wiggle output
        if towig:
            BamBins.BinsToWig(outputdir,treatments[histmark],treatBins,Binsize,chromsplit)

        # call peak
        hmmStates,lams= BamBins.callPeakFromBins(treatBins)

        # train HMM and decode HMM. hmmStates will be updated.
        m=BamBins.trainHMM(hmmStates['chr1'])
        BamBins.decodeHMM(m,hmmStates)

        # parse HMM and generate peak files. *.bam -> *.peaks
        Peaks[histmark] = BamBins.HMMToPeaks(hmmStates,treatBins,Binsize)
        print >>sys.stderr, "Write Peaks to file:",treatments[histmark].split(".")[0]+".peaks"
        fh=open(treatments[histmark].split(".")[0]+".peaks",'w')
        for tpeak in Peaks[histmark]:
            print >>fh, tpeak
        fh.close()
        print >>sys.stderr, "Write Peaks finished."
        print >>sys.stderr 
    

    


