#!/usr/bin/python
#Last-modified: 10 Jan 2013 11:14:00 AM

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
import string
from wLib.wBed import IO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------
gpos=0
gchrom=""
##
# @brief Put base into genome
#
# @param fh file handel to write the fasta sequence
# @param chrom current genome
# @param base base will be write in
#
# @return None
def putbase(fh,chrom,base=None):
    ''' Put base into fasta file.'''
    global gpos
    global gchrom
    # write header
    if gchrom!=chrom:
        if gchrom!="":
            fh.write("\n")
        fh.write(">"+chrom+"\n")
        gchrom=chrom
        gpos=0
    if base and base!='*':
        fh.write(base)
        gpos+=1
        if gpos%80==0:
            fh.write("\n")
    return None

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" *.rev prefix"
    else:
        nCount=0
        curchrom=''
        fafh=open(sys.argv[2]+".fa",'w')
        revfh=open(sys.argv[2]+".rev",'w')
        for line in IO.ColumnReader(sys.argv[1]):
            if line[0]!=curchrom: # print the header of fasta file
                curchrom=line[0]
                curpos=int(line[1])
                #print curchrom,curpos
                putbase(fafh,curchrom)
            ref=line[2] # ref base
            pos=int(line[1]) # position
            depth=int(line[3])
            basestr=line[4].upper()
            bases={'A':0,'C':0,'G':0,'T':0,'N':0,'*':0}
            i=0
            lstr=len(basestr)
            deletion=''
            insertion=''
            annotation={}
            while(i<lstr):
                if basestr[i] == '^': # head of read
                    i+=1
                elif basestr[i]== '.' or basestr[i] == ',': # match in either strand
                    bases[ref]+=1
                elif basestr[i] in "ATCGN*": # mismatch
                    bases[basestr[i]]+=1
                elif basestr[i]=='$': # end of read
                    pass
                elif basestr[i]=='-': # deletion
                    i+=1
                    deletion=int(basestr[i])
                    i+=deletion
                    annotation.setdefault(basestr[(i-deletion-1):(i+1)],0)
                    annotation[basestr[(i-deletion-1):(i+1)]]+=1
                elif basestr[i]=='+': # insertion
                    i+=1
                    insertion=int(basestr[i])
                    i+=insertion
                    annotation.setdefault(basestr[(i-insertion-1):(i+1)],0)
                    annotation[basestr[(i-insertion-1):(i+1)]]+=1
                else: # not run
                    print basestr[i]
                i+=1
            
            # find most possible base
            maxbase=ref
            for base in bases: # max base
                if bases[base]>bases[maxbase]:
                    maxbase=base
            if maxbase!=ref: # != ref
                annotation[ref+">"+maxbase]=bases[maxbase]
            
            # put 'N' in 0 depth regions
            while curpos<pos:
                nCount+=1
                curpos+=1
                putbase(fafh,curchrom,'N')
            # put base into genome
            curpos+=1
            putbase(fafh,curchrom,maxbase)
            # annotations
            lanno=''
            for i in annotation:
                lanno+=i+":"+str(annotation[i])+";"
                if "+" in i and annotation[i]>bases[maxbase]: # if most cases have insertions
                    inserts=i[1:]
                    for base in inserts:
                        putbase(fafh,curchrom,base)
            revfh.write("\t".join(line[0:4])+"\t"+maxbase+"\t"+str(bases[maxbase])+"\t"+lanno.rstrip(";")+"\n")
        fafh.close()
        revfh.close()
        print nCount

