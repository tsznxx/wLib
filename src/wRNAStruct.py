#!/usr/bin/python
#Last-modified: 25 Feb 2014 04:04:13 PM

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
# @contact: yfwang0405@gmail.com

# ------------------------------------
# python modules
# ------------------------------------

# python packages
import os
import re
import sys
import copy
import time
import fisher
import numpy
import bisect
import random
import tempfile
from itertools import izip
from subprocess import Popen, PIPE
#from scipy.cluster.vq import kmeans2
from sklearn.cluster import KMeans

# custom packages
import wRNA

# ------------------------------------
# constants
# ------------------------------------

debug = False
rctable = {'A':'U','U':'A','C':'G','G':'C'}

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

class Fasta(object):
    '''
    Fasta format
        >YIL140W
        AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
    '''
    def __init__(self,name,seq):
        ''' Initiation. '''
        self.name = name
        self.seq = seq.upper()
    def __len__(self):
        ''' length of sequence. '''
        return len(self.seq)
    def __str__(self):
        ''' Fasta in string. '''
        return ">{0}\n{1}".format(self.name,self.seq)
    def rc(self):
        '''
        Do reverse complement of self.seq. No returns. 
        '''
        self.seq = ''.join([rctable(c) for c in self.seq])

class PARS(Fasta):
    '''
    PARS format for storing sequencing depth.
    Example file:
        >tE-UUC-I
        UCCGAUAUAGUGUAACGGCUAUCACAUCACGCUUUCACCGUGGAGACCGGGGUUCGACUCCCCGUAUCGGAG
        27023.0;99.0;16.0;2.0;2.0;56.0;25.0;58.0;27.0;28.0;41.0;51.0;138.0;220.0;37.0;50.0;52.0;12.0;8.0;43.0;23.0;47.0;74.0;34.0;19.0;31.0;17.0;27.0;33.0;36.0;13.0;7.0;19.0;371.0;12.0;146.0;45.0;22.0;31.0;22.0;34.0;46.0;34.0;16.0;4.0;1.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0
        3187.0;13.0;27.0;9.0;0.0;18.0;17.0;15.0;13.0;15.0;16.0;13.0;43.0;29.0;17.0;34.0;21.0;5.0;5.0;12.0;5.0;9.0;25.0;7.0;9.0;13.0;7.0;533.0;55.0;22.0;11.0;4.0;5.0;79.0;2.0;26.0;17.0;18.0;19.0;20.0;15.0;49.0;33.0;20.0;5.0;1.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0
    '''
    def __init__(self, name, seq, loops = [], stems = []):
        ''' Initiation. '''
        super(PARS,self).__init__(name,seq)
        if len(loops) == 0:
            self.loops = numpy.zeros(len(seq),dtype=float)
        else:
            self.loops = numpy.array(loops,dtype=float)
        if len(stems) == 0:
            self.stems = numpy.zeros(len(seq),dtype=float)
        else:
            self.stems = numpy.array(stems,dtype=float)
    def __str__(self):
        ''' string of PARS format. '''
        return ">{0}\n{1}\n{2}\n{3}".format(self.name, self.seq, ';'.join([str(round(i,2)) for i in self.loops]), ";".join([str(round(i,2)) for i in self.stems]))
    def normalize(self, sratio = 1, vratio = 1, ntrim = 5):
        ''' Normalization. '''
        self.loops *= sratio
        self.stems *= vratio
        if ntrim == 0:
            return
        self.loops[0:ntrim] = 0
        self.stems[0:ntrim] = 0
        self.loops[-ntrim:] = 0
        self.stems[-ntrim:] = 0
        return
    def rc(self):
        ''' Reverse complement. '''
        super(PARS,self).rc()
        self.loops = self.loops[::-1]
        self.stems = self.stems[::-1]
        return
    def toFastC(self, sthreshold = 2, vthreshold = -2):
        ''' Calculate constraints from PARS depth. '''
        pars = numpy.log2((self.loops+1)/(self.stems+1))
        constraints = numpy.repeat('.', len(self.seq))
        constraints[numpy.where(pars >  sthreshold)] = 'x'
        constraints[numpy.where(pars <  vthreshold)] = '|'
        return FastC(self.name, self.seq, constraints.tostring())
    def toFastCByFisher(self, Ssum, Vsum, threshold = 0.05):
        ''' Fisher exact test to calculate the pvalue and then calculate fdr. '''
        SVsum = (Ssum+Vsum)/2.
        pvalues = numpy.apply_along_axis(lambda x: fisher.pvalue(x[0],Ssum-x[0],x[1],Vsum-x[1]),1,zip(self.loops,self.stems))

class FastC(Fasta):
    '''
    FastC format for RNA structure prediction.
    Example file:
        >YIL140W
        AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
        .|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........
    Note:
        '.' means free to fold
        'x' means loops
        '|' means stems

    '''
    def __init__(self, name, seq, constraints=None):
        super(FastC,self).__init__(name,seq)
        self.constraints = constraints
        if constraints is None:
            self.constraints = '.'*len(self.seq)
    def __str__(self):
        ''' FastC string. '''
        return ">{0}\n{1}\n{2}".format(self.name,self.seq,self.constraints)

class FastS(Fasta):
    ''' 
    RNA structure.
    Example: Multiple structures are supported.
        >YIL140W
        AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
        .......((((((((............)))).))))..............  (-6.0)
        .......(((((((..............))).))))..............  (-4.0)
    '''
    def __init__(self, name, seq, structures, scores=None):
        '''
        Initiation.

        Parameters:
            name: string
                Name of sequence.
            seq: string
                RNA sequences.
            structures: string or list of strings
                RNA structure in Dot-Bracket format.
            scores: None, float or list of float
                A score of each structure. Could be the MFE score, the percentage or any other scores.
        Returns:
            fs: FastS object
                FastS object
        '''
        super(FastS,self).__init__(name,seq)
        if isinstance(structures,basestring):
            self.structures = [structures]
        else:
            self.structures = [st for st in structures]
        if scores is None:
            self.scores = [1./len(self.structures)] * len(self.structures)
        elif isinstance(scores,int) or isinstance(scores,float):
            self.scores = [float(scores)]
        else:
            self.scores = [float(sc) for sc in scores]
    def __str__(self):
        ''' FastS string. '''
        lstr = ">{0}\n{1}".format(self.name,self.seq)
        for st, sc in izip(self.structures, self.scores):
            lstr += "\n{0}\t({1})".format(st, round(sc,3))
        return lstr

class IO(object):
    ''' Process raw data to FastC format. '''
    def readDataToPARS(fafile,S,V,maxlen=300,ntrim=5,threshold = 0.05, prefix=None):
        ''' read fasta, S1 and V1 files and normalize by read number. '''
        # print parameters
        print >>sys.stderr, "Parameters:"
        print >>sys.stderr, "\tFasta file:", fafile
        print >>sys.stderr, "\tS1 file:", S
        print >>sys.stderr, "\tV1 file:", V
        print >>sys.stderr, "\tMaximum sequence length:", maxlen > 0 and maxlen or "no-limit"
        print >>sys.stderr, "\tTrim {0} nt at both ends".format(ntrim)
        print >>sys.stderr, "\tThreshold for constraints:", threshold
        if prefix is None:
            prefix = os.path.splitext(os.path.basename(fafile))[0]
        print >> sys.stderr, "\tOutput file prefix:", prefix
        print >>sys.stderr

        # read S1 file
        print >>sys.stderr, "Reading {0}".format(S)
        starttime=time.clock()
        Ssum = 0
        S1 = {}
        i0 = 0
        with open(S) as fh:
            for s in fh:
                i0 += 1
                if i0 % 100 == 0:
                    print >> sys.stderr, "Prcessed ",i0, "lines.     \r",
                items = s.split()
                S1[items[0]] = [float(i) for i in items[2].split(";")]
                Ssum += sum(S1[items[0]][ntrim:-ntrim])
        print >> sys.stderr, "Prcessed ",i0, "lines.     "
        print >>sys.stderr, "Time used: %-3.2f secs.\n" % (time.clock()-starttime)

        # Read V1 file
        print >>sys.stderr, "Reading {0}".format(V)
        starttime=time.clock()
        Vsum = 0
        V1 = {}
        i0 = 0
        with open(V) as fh:
            for v in fh:
                i0 += 1
                if i0 % 100 == 0:
                    print >> sys.stderr, "Prcessed ",i0, "lines.     \r",
                items = v.split()
                V1[items[0]] = [float(i) for i in items[2].split(";")]
                Vsum += sum(V1[items[0]][ntrim:-ntrim])
        print >> sys.stderr, "Prcessed ",i0, "lines.     "
        print >>sys.stderr, "Time used: %-3.2f secs.\n" % (time.clock()-starttime)
        
        # Calculating normalization factors
        print >>sys.stderr, "Calculating normalization factors:"
        print >>sys.stderr, "Number of sequences:",min(len(V1),len(S1))
        print >> sys.stderr, "sum(S) = {0}, sum(V) = {1}".format(Ssum, Vsum)
        sk = (Ssum + Vsum)/2.0/Ssum
        vk = (Ssum + Vsum)/2.0/Vsum
        print >> sys.stderr, "Normalization ratio: S = %-3.2f, V = %-3.2f\n" % (sk, vk)

        # Read fasta file and generate PARS
        print >> sys.stderr, "Readling fasta file and generating PARS list ..."
        starttime=time.clock()
        PARSlst = []
        FastClst = []
        i0 = 0
        with open(fafile) as fh, open(prefix+".fp", 'w') as fhp, open(prefix+".fc", 'w') as fhc:
            for fa in IO.FastaReader(fafile):
                i0 += 1
                if i0 % 100 == 0:
                    print >> sys.stderr, "Prcessed ",i0, "items.     \r",
                if maxlen == 0 or len(fa.seq) <= maxlen:
                    tfp = PARS(fa.name,fa.seq,S1[fa.name],V1[fa.name])
                    # p value and fdr
                    pv = [fisher.pvalue(*x) for x in zip(tfp.loops,Ssum-tfp.loops,tfp.stems,Vsum-tfp.stems)]
                    fdr_S = Algorithm.fdr_bh([i.right_tail for i in pv])
                    fdr_V = Algorithm.fdr_bh([i.left_tail  for i in pv])
                    # calculate constraints
                    cons = numpy.repeat(".",len(tfp))
                    for i in range(ntrim,len(cons)-ntrim): # flanking ntrim bases are not considered.
                        if fdr_S[i] < threshold:
                            cons[i] = "x"
                        elif fdr_V[i] < threshold:
                            cons[i] = '|'
                    # Normalization
                    tfp.normalize(sk, vk, ntrim)
                    PARSlst.append(tfp)
                    print >> fhp, tfp

                    # print constraints file
                    tfc = FastC(tfp.name,tfp.seq,cons.tostring())
                    FastClst.append(tfc)
                    print >> fhc, tfc

        print >> sys.stderr, "Prcessed ",i0, "items.     "
        print >> sys.stderr, "{0} sequences passed the maxlen filter.".format(len(PARSlst))
        print >>sys.stderr, "Time used: %-3.2f secs.\n" % (time.clock()-starttime)

        print >> sys.stderr, "Done."
        return (Ssum, Vsum, PARSlst, FastClst)
    readDataToPARS=staticmethod(readDataToPARS)
    def FastaReader(infile):
        '''Read sequence files.'''
        # Read lines
        with open(infile) as fh:
            line = fh.next()
            if line[0] != ">":
                raise ValueError("Records in Fasta files should start with '>' character")
            line = line.lstrip('>').split()
            name = line[0]
            seq = ''
            while True:
                try:
                    line = fh.next()
                except:
                    if seq != '':
                        yield Fasta(name, seq)
                    raise StopIteration
                if line[0] != ">":
                    seq += line.rstrip()
                else:
                    yield Fasta(name, seq)
                    line = line.lstrip('>').split()
                    name = line[0]
                    seq = ''
    FastaReader=staticmethod(FastaReader)
    def PARSReader(infile):
        ''' Read PARS file (sequences, S1 depth and V2 depth) and calculate the constraints with threshold. '''
        with open(infile) as fh:
            for line in fh:
                name = line.lstrip(">").rstrip()
                seq  = fh.next().rstrip()
                S    = numpy.array([float(i) for i in fh.next().split(";")])
                V    = numpy.array([float(i) for i in fh.next().split(";")])
                yield PARS(name, seq, S, V)
        raise StopIteration
    PARSReader=staticmethod(PARSReader)
    def FastCReader(infile):
        ''' Read FastC file. '''
        with open(infile) as fh:
            for line in fh:
                name = line.lstrip(">").rstrip()
                seq  = fh.next().rstrip()
                constraints = fh.next().rstrip()
                yield FastC(name, seq, constraints)
    FastCReader=staticmethod(FastCReader)
    def FastSReader(infile):
        ''' Read FastS file. '''
        with open(infile) as fh:
            name = ''
            seq = ''
            structures = []
            scores = []
            for line in fh:
                if line.startswith('>'):
                    if name != '':
                        yield FastS(name, seq, structures, scores)
                    name = line.lstrip('>').rstrip()
                    seq  = fh.next().rstrip()
                    structures = []
                    scores = []
                else:
                    lstr = line.split()
                    structures.append(lstr[0])
                    try:
                        scores.append(float(lstr[1].rstrip(")").lstrip("(")))
                    except:
                        scores.append(1.0)
            if name != '':
                yield FastS(name, seq, structures, scores)
        raise StopIteration
    FastSReader=staticmethod(FastSReader)

class Predictor(object):
    '''
    RNA structure predictor.
    Usage:
        structure, score = Predictor.RNAfold(tFastC, Temperature, threshold)
        print score, structure
        
    '''
    home = os.path.expanduser('~')
    def FoldMerge(tfastc, T = 37, predictors=['RNAfold','Fold','pknots','UNAFold','sfold']):
        ''' Merge the folded structures from multiple software. '''
        predictor_dict = {'RNAfold': Predictor.RNAfold, 'Fold': Predictor.Fold, 'pknots':  Predictor.pknots, 'mfold':  Predictor.mfold, 'UNAFold':  Predictor.UNAFold, 'sfold':Predictor.sfold,'pknotsRG':Predictor.pknotsRG, 'ipknot':Predictor.ipknot}
        structures = []
        scores = []
        for predictor in predictors:
            if predictor_dict.has_key(predictor):
                tfasts = predictor_dict[predictor](tfastc)
                structures.extend(tfasts.structures)
                scores.extend(tfasts.scores)
        # unique structure
        struct_dict= {}
        for st, sc in izip(structures,scores):
            if struct_dict.has_key(st):
                struct_dict[st] = min(sc, struct_dict[st]) # energy the lower the better
            else:
                struct_dict[st] = sc
        return FastS(tfastc.name, tfastc.seq, struct_dict.keys(), struct_dict.values())
    FoldMerge=staticmethod(FoldMerge)
    def RNAfold(tfastc, T = 37, threshold = 0):
        ''' Call RNAfold to predict the structure. '''
        # Run RNAfold with contraints
        p = Popen(['RNAfold','-C', '-T', str(T)], stdin = PIPE, stdout = PIPE, stderr = PIPE)
        p.stdin.write(str(tfastc))
        pstdout, pstderr = p.communicate()
        
        # Check the threshold here
        # to do 

        # Extract structure and score
        pre = re.compile(r'(\S+)\s+\((\s*-*\d+\.\d+)\)')
        s = pstdout.rstrip().split('\n')[-1]
        structure , score = pre.search(s).groups()
        if os.path.isfile(tfastc.name+"_ss.ps"):
            os.remove(tfastc.name+"_ss.ps")
        return FastS(tfastc.name, tfastc.seq, [structure], [float(score)])
    RNAfold=staticmethod(RNAfold)
    def Fold(tfastc, T = 37, threshold = 0):
        '''
        Call Fold (RNA Structure package) to predict structures.
        Note: Fold won't give structures if it cannot satisfy all the constraints. [Fix this in the future].
        '''
        # check DATAPATH environment variable
        if not os.environ.has_key('DATAPATH'):
            raise KeyError("Fold (RNA Structure) Error: Please set environment variable $DATAPATH to the location of the data_tables.")
        T = Utils.TempConverter(T,'C','K') # from C to K
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa')
        (fdct, ctname) = tempfile.mkstemp(suffix='.ct')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.CON')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        with os.fdopen(fdcn, 'w') as fh:
            fh.write(Utils.PARSToFold(tfastc.constraints))
        os.close(fdct)
        structures = []
        scores = []
        try:
            # Run Fold
            # Fold 
            p = Popen(['Fold', faname, ctname, '-T', str(T), '-C', cnname], stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            if pstderr:
                raise ValueError("Fold (RNA Structure package) run error: {0}".format(pstderr))
            # Parse the *.ct file
            tfs = Utils.ct2dot(ctname)
        finally:
            os.remove(faname)
            os.remove(ctname)
            os.remove(cnname)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by Fold."
        return tfs # FastS(tfastc.name, tfastc.seq, structures, scores)
    Fold=staticmethod(Fold)
    def pknots(tfastc, T = 37, threshold = 0):
        '''
        Fold by pknots. 
        pknots doesn't have temperature parameter.
        pknots has -k parameter for pseudoknots folding.
        '''
        if T != 37:
            print >> sys.stderr, "WARNING: pknots doesn't accept temperature != 37 !"
            return FastS(tfastc.name, tfastc.seq, [], [])
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa')
        (fdct, ctname) = tempfile.mkstemp(suffix='.ct')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}".format(tfastc.name,tfastc.seq))
        os.close(fdct)
        structures = []
        scores = []
        try:
            # Run pknots
            p = Popen(['pknots', '-k', faname, ctname], stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            if pstderr:
                raise ValueError("ERROR: pknots run error: {0}".format(pstderr))
            # Parse the *.ct file
            #structures, scores = Utils.ct2dot(ctname)
            structures, scores = Utils.pknots2dot(ctname)
        finally:
            os.remove(faname)
            os.remove(ctname)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by pknots."
        return FastS(tfastc.name, tfastc.seq, structures, scores)
    pknots=staticmethod(pknots)
    def pknotsRG(tfastc, T = 37, threshold = 0):
        ''' Fold by pknotsRG. '''
        if T != 37:
            print >> sys.stderr, "WARNING: pknotsRG doesn't accept temperature != 37 !"
            return FastS(tfastc.name, tfastc.seq, [], [])
        structures = []
        scores = []
        try:
            p = Popen(['pknotsRG', '-s'], stdin = PIPE, stdout = PIPE, stderr = PIPE)
            p.stdin.write(">{0}\n{1}".format(tfastc.name, tfastc.seq))
            pstdout, pstderr = p.communicate()
        except OSError:
            raise OSError('pknotsRG is not installed.')
        if pstderr:
            raise IOError('ERROR: pknotsRG run error: '+pstderr)
        # parse the output
        for line in pstdout.split('\n')[6:-1]: # first 5 lines are descriptions
            structure, score = line.split()
            score = float(score[1:-2]) ## (-5.80)
            if score < 0 and structure not in structures:
                structures.append(structure)
                scores.append(score) 
        return FastS(tfastc.name, tfastc.seq, structures, scores)
    pknotsRG=staticmethod(pknotsRG)
    def mfold(tfastc, T = 37, threshold = 0):
        '''
        Fold by mfold.
        mfold SEQ=test.fa AUX=test.con T=37 RUN_TYPE=html
        '''
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.con')
        with os.fdopen(fdfa, 'w') as fh: # with statement will close the fh and fd in the end
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        with os.fdopen(fdcn, 'w') as fh:
            fh.write(Utils.PARSToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # mfold 
            p = Popen(['mfold', 'SEQ='+faname, 'AUX='+cnname, 'T={0}'.format(T), 'RUN_TYPE=html'], stdin = None, stdout = PIPE, stderr = PIPE) 
            pstdout, pstderr = p.communicate()
            if pstderr:
                raise ValueError("ERROR: mfold run error: {0}".format(pstderr))
        finally:
            os.remove(faname)
            os.remove(cnname)
            # Parse the *ct file
            idx = 1
            bname = os.path.basename(faname)
            fn = "{0}_{1}.ct".format(bname,idx)
            while os.path.isfile(fn):
                tfs = Utils.ct2dot(fn)
                structures.extend(tfs.structures)
                scores.extend(tfs.scores)
                idx += 1
                fn = "{0}_{1}.ct".format(bname,idx)
            for fn in os.listdir("."):
                if fn.startswith(bname):
                    os.remove(fn)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by mfold."
        if os.path.isfile('date.test'):
            os.remove('date.test')
        return FastS(tfastc.name, tfastc.seq, structures, scores) 
    mfold=staticmethod(mfold)
    def UNAFold(tfastc, T = 37, threshold = 0):
        ''' Fold by UNAFold. '''
        # Create temp file for fa and constraints
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa',dir=".")
        (fdcn, cnname) = tempfile.mkstemp(suffix='.aux',dir=".")
        faname = os.path.basename(faname)
        cnname = os.path.basename(cnname)
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        with open(cnname, 'w') as fh:
            fh.write(Utils.PARSToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # mfold 
            p = Popen(['UNAFold.pl', '-t', str(T), '-c', cnname, '--run-type=html', faname], stdin = None, stdout = PIPE, stderr = PIPE) 
            pstdout, pstderr = p.communicate()
            if pstderr:
                print >> sys.stderr, "ERROR: UNAFold run error: {0}".format(pstderr)
        finally:
            os.remove(faname)
            os.remove(cnname)
            # Parse the *ct file
            bname = os.path.basename(faname)
            for fn in os.listdir("."):
                if fn.startswith(bname):
                    if fn.endswith(".ct"):
                        tfs = Utils.ct2dot(fn)
                        for st,sc in izip(tfs.structures,tfs.scores):
                            if st not in structures:
                                structures.append(st)
                                scores.append(sc)
                    os.remove(fn)
        if len(structures) == 0:
            print >> sys.stderr, "No structures predicted by UNAFold."
        return FastS(tfastc.name, tfastc.seq, structures, scores) 
    UNAFold=staticmethod(UNAFold)
    def ipknot(tfastc, T = 37, threshold = 0):
        ''' ipknot fold. '''
        if T != 37:
            print >> sys.stderr, "WARNING: pknotsRG doesn't accept temperature != 37 !"
            return FastS(tfastc.name, tfastc.seq, [], [])
        (fdfa, faname) = tempfile(suffix='.fa')
        with os.fdopen(fdfa) as fh:
            fh.write(">{0}\n{1}".format(tfastc.name, tfastc.seq))
        try:
            p = Popen(['ipknot', '-r', '10', faname], stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
        except OSError:
            raise OSError('ERROR: ipknot is not installed.')
        finally:
            os.remove(faname)
        # Parse the result
        if pstderr:
            raise IOError('ERROR: ipknot run error: '+pstderr)
        structure =  pstdout.split('\n')[-2].split()
        if structure.strip(".") != "":
            return FastS(tfastc.name, tfastc.seq, [structure], [0.0])
        return FastS(tfastc.name, tfastc.seq, [], [])
    ipknot=staticmethod(ipknot)
    def sfold(tfastc, T = 37, threshold = 0):
        ''' Fold by UNAFold. '''
        # Create temp file for fa and constraints
        prefix = tempfile.mkdtemp(dir=".")
        (fdfa, faname) = tempfile.mkstemp(suffix='.fa',dir='.')
        (fdcn, cnname) = tempfile.mkstemp(suffix='.aux',dir='.')
        with os.fdopen(fdfa, 'w') as fh:
            fh.write(">{0}\n{1}\n".format(tfastc.name,tfastc.seq))
        # constraints
        cns = False
        if tfastc.constraints.find('|') != -1 or tfastc.constraints.find('x') != -1:
            cns = True
            with os.fdopen(fdcn, 'w') as fh:
                fh.write(Utils.PARSToMfold(tfastc.constraints))
        structures = []
        scores = []
        try:
            # sfold
            if cns:
                p = Popen([Predictor.home+'/software/RNAStructure/sfold-2.2/bin/sfold', '-f', cnname, '-o', prefix, faname], stdin = None, stdout = PIPE, stderr = PIPE)
            else:
                p = Popen([Predictor.home+'/software/RNAStructure/sfold-2.2/bin/sfold', '-o', prefix, faname], stdin = None, stdout = PIPE, stderr = PIPE)
            pstdout, pstderr = p.communicate()
            if pstderr:
                print >> sys.stderr, "ERROR: sfold run error: {0}".format(pstderr)
            if os.path.isdir(prefix+"/clusters"): # clusters exists
                for f in os.listdir(prefix+"/clusters"):
                    if f.endswith(".ct"):
                        tfs = Utils.ct2dot(prefix+"/clusters/"+f)
                        for st,sc in izip(tfs.structures,tfs.scores):
                            structures.append(st)
                            scores.append(sc)
        finally:
            os.remove(faname)
            os.remove(cnname)
            os.system('rm -rf '+prefix)
        return FastS(tfastc.name,tfastc.seq, structures, scores)
    sfold=staticmethod(sfold)
    def snoGPS(tfa, T = 37, threshold =0):
        '''
        Predict H/ACA box snoRNAs.
        '''
        pass
        return tfs
    snoGPS=staticmethod(snoGPS)

class Algorithm(object):
    ''' Algorithm for RNA structure prediction. '''
    def Greedy(tfastpars, tfasts, sthreshold = 2, vthreshold = -2):
        '''
        Find the structures at match the experiment data.
        Example:
            loops       = 1  3  3  0  1  6  6  0  4  4  0  3  3  0  1  6  6  0  1  
            stems       = 0  1  1  3  3  1  0  3  0  1  3  0  1  6  6  0  1  1  0 
            cons        = .  x  x  |  |  x  x  |  x  x  |  x  x  |  |  x  x  .  .
            struct      = .  .  .  (  (  (  .  (  .  .  )  .  )  )  )  .  .  .  . 
            loop score       3  3           6     4  4     3           6  6  
            stem score             3  3        3        3        6  6 
        The score for this strtucture will be sum of loop and stem depth that satisfy the constraints.
        Then we remove the depth compatible to the best structure by taking their median depth out.
        loops = loops - 4 and stems = stems - 3
            loops       = 1  0  0  0  1  6  2  0  0  0  0  0  3  0  1  2  2  0  1  
            stems       = 0  1  1  0  0  1  0  0  0  1  0  0  1  3  3  0  1  1  0 
            cons        = .  .  .  |  |  x  x  |  x  x  |  x  x  |  |  x  x  .  .
            struct      = .  (  (  .  (  .  .  (  .  .  )  .  )  )  )  .  .  .  . 
            loop score       3  3           6     5  4     3           6  6  
            stem score             3  3        3        3        6  6 

        '''
        # single structure
        if len(tfasts.structures) ==1:
            return copy.deepcopy(tfasts)
        # multiple structures
        # data preparation
        loops = copy.deepcopy(tfastpars.loops)
        stems = copy.deepcopy(tfastpars.stems)
        tfastc = tfastpars.toFastC(sthreshold, vthreshold)
        cons = numpy.array(list(tfastc.constraints))
        structures = []
        for st in tfasts.structures:
            structures.append(numpy.array(list(st)))
        scores = copy.deepcopy(tfasts.scores)
        # start iteration
        matchstructures = []
        matchscores = []
        medianscore = 1000000
        while medianscore > 0 and len(structures) > 0:

            # find best match
            bestscore = 0
            bestidx = 0
            bestenergy = scores[0]
            for i,st in enumerate(structures):
                lscore = numpy.median(loops[numpy.where(st == '.')])
                sscore = numpy.median(stems[numpy.where(st != '.')])
                score = max(lscore, sscore)
                if score > bestscore or (score == bestscore and scores[i] < bestenergy): # best score or (equal score and best energy)
                    bestidx = i
                    bestscore = score
                    bestenergy = scores[i]
            medianscore = bestscore
            print bestscore
            # remove the compatible depth [ to be modified in the future.]
            loops[numpy.where(structures[i] == '.')] -= numpy.median(loops[numpy.where(structures[i] == '.')])
            loops[loops<0] = 0
            stems[numpy.where(structures[i] != '.')] -= numpy.median(stems[numpy.where(structures[i] != '.')])
            stems[stems<0] = 0
            # record the best structure in the matchstructures
            matchstructures.append(structures[i].tostring())
            matchscores.append(scores[i])
            # remove the best structure from structures.
            del structures[i]
            del scores[i]
        return FastS(tfasts.name,tfasts.seq, matchstructures, matchscores)
    Greedy=staticmethod(Greedy)
    def simulation(fs,numR,nS=1,nV=1):
        '''
        Simulate reads for structures with a given percentages.
            fs: FastS object contains multiple structures and an arbitrary proprotion.
            numR: number of RNA molecular
            nS, nV: proportion of reads for S and V.
        Note:
            numR is attributed to loops and stems by nS:nV.
            nS and nV are attributed to structures by scores of structures.
        '''        
        # initiation
        L = len(fs)
        M = len(fs.structures)
        S = numpy.zeros(L)
        V = numpy.zeros(L)
        # ranges of each structure
        scores = numpy.array(fs.scores,dtype=float)
        # Number of RNA molecular for S and V
        nS = round(1.0*numR*nS/(nS+nV))
        nV = numR - nS
       
        # z matrix
        z = numpy.zeros((L,M))
        for k in range(L):
            for j in range(M):
                if fs.structures[j][k] != ".": # k is stem in j
                    z[k,j] = 1
                    
        # generate random reads
        numpy.random.seed(int(time.time()*100)%100)
        for j in xrange(M):
            # nuber of RNA molecular for each structure
            # loops
            idx = numpy.nonzero(1-z[:,j])[0]
            l = len(idx)
            nSj = round(nS*scores[j]/sum(scores))*l
            for i in numpy.random.randint(l,size=nSj):
                S[idx[i]] += 1
            # stems
            idx = numpy.nonzero(z[:,j])[0]
            l = len(idx)
            nVj = round(nV*scores[j]/sum(scores))*l
            for i in numpy.random.randint(l,size=nVj):
                V[idx[i]] += 1
        return PARS(fs.name,fs.seq,S,V)
    simulation=staticmethod(simulation)
    def EM(fp,fs,centroids=None,threshold = 1e-6, maxiter=100):
        ''' 
        EM algorithm to calculate a percentage for each structure centroid. 
        '''
        def normalize(x,p):
            s = sum(x*p)
            if s == 0:
                return x
            else:
                return x*p/s
        
        if centroids is None:
            centroids = range(len(fs.structures))
        M = len(centroids) # Number of clusters
        if M ==0:
            return [0],[0,0,0,0]
        L = len(fp.loops)  # Length of RNA
        
        S = numpy.zeros(L)
        V = numpy.zeros(L)
        S += fp.loops
        V += fp.stems
       
        # Calculate Zkj
        z = numpy.zeros((L,M))
        for k in range(L):
            for j in range(M):
                if fs.structures[centroids[j]][k] != ".": # k is stem in j
                    z[k,j] = 1
        # Initial percentage
        pi = numpy.zeros(M)+1.0/M
        zV = numpy.apply_along_axis(sum,0,z) # number of stems in each structure
        zS = L-zV # number of loops in each structure
        
        # iteration
        cnt = 0
        smS = 0
        smV = 0
        while (True):
            # E step: E[P(X,Z|theta)]
            # normalize the possibilities by row
            ES = numpy.apply_along_axis(lambda x:normalize(x,pi),1,1-z)
            EV = numpy.apply_along_axis(lambda x:normalize(x,pi),1,z)
            # expected reads for each structure
            muS = numpy.zeros(M)
            muV = numpy.zeros(M)
            for k in xrange(L):
                for j in xrange(M):
                    muS[j] += S[k]*ES[k,j]
                    muV[j] += V[k]*EV[k,j]
            smS = sum(muS)
            smV = sum(muV)
            mu = muS/zS + muV/zV
            # M step: argmax(E[P(X,Z|theta)])
            mu /= sum(mu)
            
            # check if converge
            if cnt > maxiter or numpy.sum(numpy.abs(mu-pi)) < threshold:
                pi = mu
                break
            pi = mu
            cnt += 1
        print >> sys.stderr, "EM iterates {0} times.".format(cnt)
        return pi, [smS,sum(S),smV,sum(V)]
    EM=staticmethod(EM)
    def fitness(fp,fs,pi=[1]):
        '''
        Calculate fitness, the percentage of compatible reads, given all reads, structures and percentages: 
        Input: 
            fp: normalized PARS object. Flanking bases and outliers are trimmed.
            fs: FastS object with known or predicted structures.
            pi: percentage of each structures. sum(pi)=1.
        Output:
            smS: number of compatible S reads.
            S:   total number of S reads.
            smV: number of compatible V reads.
            V:   total number of V reads.
        '''
        def normalize(x,p):
            ''' Normalize conditional probabilities. '''
            s = sum(x*p)
            if s == 0:
                return x
            else:
                return x*p/s
        M = len(fs.structures)
        L = len(fp)
        if M == 0:
            return [0,0,0,0]
        S = numpy.zeros(L)
        V = numpy.zeros(L)
        S += fp.loops
        V += fp.stems
        # Calculate Zij
        z = numpy.zeros((L,M))
        for i in range(L):
            for j in range(M):
                if fs.structures[j][i] != ".": # k is stem in j
                    z[i,j] = 1
        # Expected depth
        ES = numpy.apply_along_axis(lambda x:normalize(x,pi),1,1-z)
        EV = numpy.apply_along_axis(lambda x:normalize(x,pi),1,z)
        # expected reads for each structure
        muS = numpy.zeros(M)
        muV = numpy.zeros(M)
        for i in xrange(L):
            for j in xrange(M):
                muS[j] += S[i]*ES[i,j]
                muV[j] += V[i]*EV[i,j]
        smS = sum(muS)
        smV = sum(muV)
        return smS,sum(S),smV,sum(V)
    fitness=staticmethod(fitness)
    def spectralClustering(data):
        '''
        Do spectral clustering on symmetric distance matrix, data, into K groups, in which K is determined by eigengap heuristic method [1].
        Usage:
            centroids, clusters = spectralClustering(data,k)
        Example:
            data = [[ 0.  4.  4.  2.  2.]
                    [ 4.  0.  4.  4.  4.]
                    [ 4.  4.  0.  4.  4.]
                    [ 2.  4.  4.  0.  4.]
                    [ 2.  4.  4.  4.  0.]]
            k = 2
            centroid_idx, clusters = spectralClustering(data)
            Output:
                > print centroids
                [0,1]
                > print clusters
                [array([0, 3, 4]), array([1, 2])]
        References:
            [1] von Luxburg, U. (2007). "A tutorial on spectral clustering." Statistics and Computing 17(4): 395-416.
        '''
        n = data.shape[0]
        sd = numpy.std(data)
        
        # Affinity matrix
        A = numpy.zeros(data.shape,dtype=float)
        for i in range(n):
            for j in range(i):
                A[i,j] = numpy.exp(-sum((data[i,:]-data[j,:])**2/2/sd**2))
                A[j,i] = A[i,j]
                
        # Diagonal matrix
        D = numpy.zeros(data.shape)
        for i in range(n):
            D[i,i] = sum([A[i,j] for j in range(n)])
            
        # Scaled matrix L from D
        S = numpy.linalg.inv(D)**0.5
        L = numpy.dot(S,A)
        L = numpy.dot(L,S)
        # Eigen values and vectors: w, v
        w, v = numpy.linalg.eig(L)
        
        # find the best K using eigengap heuristic method 
        # von Luxburg, U. (2007). "A tutorial on spectral clustering." Statistics and Computing 17(4): 395-416.
        sw = w.argsort()
        eigengaps = [w[sw[i+1]]-w[sw[i]] for i in range(n-1)]
        k = numpy.argmax(eigengaps) + 1
        # plot(range(n-1),eigengaps,'*')
        
        # get the eigen vectors for the first smallest eigen values    
        X = v[:,sw[:k]]

        # normalize X
        Y = numpy.apply_along_axis(lambda x: x/numpy.sqrt(sum(x**2)),1,X)
        
        # do k-means clustering
        km = KMeans(n_clusters=k, init='k-means++', max_iter=100, n_init=5, verbose=False)
        rst = km.fit(Y)
        centroids = rst.cluster_centers_
        labels = rst.labels_
        
        # index of centroids
        centroid_idx = numpy.zeros(k,dtype=numpy.int32)
        for i,centroid in enumerate(centroids):
            cur = 0
            curdis = -1
            for j,point in enumerate(Y):
                dis = numpy.sum((centroid-point)**2)
                if dis < curdis or curdis == -1:
                    cur = j
                    curdis = dis
            centroid_idx[i]=cur
        return centroid_idx,labels
    spectralClustering=staticmethod(spectralClustering)    
    def kmeans(data, k, maxiter = 20 ): # obsoleted
        ''' K-means clustering. '''
        n = data.shape[0] 
        labels = numpy.repeat(0,n)
        minima = None # global minima
        for i in xrange(maxiter):
            # generate k centers randomly
            centroids = random.sample(range(n),k) # current k centers
            if debug: print "Start K-means: k =", k, ", # of samples:", n
            while True:
                # get labels for each elements    
                for j in range(n):
                    mindis = data[j,centroids[0]] # set label to the first centroid
                    labels [j] = centroids[0]
                    for i in range(k): # find the nearest centroid
                        if data[j,centroids[i]] < mindis:
                            labels[j] = centroids[i]
                            mindis = data[j,centroids[i]] # the nearest centroid
                if debug: print "centroids", centroids, "labels",labels
                # Find new centroids
                precentroids = [i for i in centroids]
                distance_square_sum = numpy.repeat(10000000,k) # distance square sum 
                # For each point, check if it is the centroid of its cluster.
                for i in range(n):
                    sumi = 0. # sum of square distance
                    cnt = 0 # count of samples
                    for j in range(n):
                        if labels[i] == labels[j]: 
                            sumi += data[i,j]**2
                            cnt  +=1
                    sumi /= cnt
                    idx = precentroids.index(labels[i]) # index of the centroid/ksum 
                    if sumi < distance_square_sum[idx]: # possible new center
                        distance_square_sum[idx] = sumi
                        centroids[idx] = i
                if debug: print "Total score:",sum(distance_square_sum)
                if numpy.array_equal(centroids,precentroids): # check if converge
                    break
            # Check if minima.
            tsum = sum(distance_square_sum)
            if minima is None or tsum < minima:
                minima = tsum
                bestcluster = (sorted(centroids),numpy.array(labels),tsum)
        return bestcluster
    kmeans=staticmethod(kmeans)
    def BIC(data, centroids, labels, seqlen): # not working for unknown reason
        ''' BIC algorithm to find the best k of K-means. '''
        n = len(labels)
        k = len(centroids)
        R = len(labels)
        M = seqlen # M should be the dimension of the points.
        Ri = [len(labels[labels == centroids[i]]) for i in range(k)]
        if R == k:
            sigma = 0
        else:
            sigma = sum([data[i,labels[i]]**2 for i in xrange(R)])/float(R-k)
        lD = -0.5*numpy.log(2*numpy.pi) - 0.5*R*M*sigma - 0.5*R + 0.5*k**2 + sum([r*numpy.log(r) for r in Ri]) - R*numpy.log(R)
        return lD - 0.5*(M+1)*k*numpy.log(R)
    BIC=staticmethod(BIC)
    def hcluster(ndata,threshold=20):
        ''' Hierarchical clustering. '''
        data = numpy.array(ndata)
        R = data.shape[0]
        clusters = {i:[i] for i in range(R)} # cluster initiation.
        for k in range(R-1):
            # Find the pair with the shortest distance
            idx = numpy.where(data == numpy.min(data[data!=0]))
            i = idx[0][0]
            j = idx[1][0]
            if debug: print i,j,data[i,j]
            if data[i,j] <= threshold:                
                # Update data
                if debug: print data[i,:],len(clusters[i]),len(clusters[j])
                data[i,:] = (data[i,:]*len(clusters[i]) + data[j,:]*len(clusters[j]))/(len(clusters[i])+len(clusters[j]))
                data[:,i] = data[i,:]
                clusters[i].extend(clusters[j])
                del clusters[j]
                data[i,i] = 0
                data[j,:] = 0
                data[:,j] = 0
                if debug: print data
            else: # No pairs have distance < threshold
                break
        return clusters.values()
    hcluster=staticmethod(hcluster)    
    def fdr_bh(pvalues):
        n = len(pvalues)
        if n == 0:
            return []
        # sort pvlaues
        args, pv = zip(*sorted(enumerate(pvalues),key= lambda x:x[1]))
        bh_values = [0] * n
        pre_v = 0
        for i, pv in enumerate(pv):
            cur_v = min(1, pv * n /(i+1)) # correction
            pre_v = max(cur_v, pre_v)
            bh_values[args[i]] = pre_v
        return bh_values
    fdr_bh=staticmethod(fdr_bh)
    def rank(x):
        ''' rank of a list. '''
        array=numpy.array(x)
        temp = array.argsort()
        ranks = numpy.arange(len(array))[temp.argsort()]
        return ranks
    rank=staticmethod(rank)
    def trimed(x, lowbound = 0.05, upbound = 0.95):
        ''' Trim the top ones and the bottom ones, and set them to bound values. '''
        r = Algorithm.rank(x) # type(r) is numpy.ndarray
        l = len(r) - 1
        ub = x[numpy.where(r == round(upbound * l))[0]]
        lb = x[numpy.where(r == round(lowbound *l))[0]]
        if type(x) == numpy.ndarray:
            nx = numpy.array(x)
            nx[nx>ub] = ub
            nx[nx<lb] = lb
        else:
            nx = [ ub if i > ub else i for i in x]
            nx = [ lb if i < lb else i for i in nx]
        return nx
    trimed=staticmethod(trimed)
    def distanceMatrix(fs,method="rbp"):
        '''
        Calculate the distance matrix, given a FastS object with multiple structures.
        Options:
            method: 'top' means calculate the topological distance using wRNA.distance from Vienna RNA package [1].
                    'rbp'  means calculate the relaxed base-pair distance [2].
        References:
            [1] Hofacker, I. L., et al. (1994). "Fast Folding and Comparison of Rna Secondary Structures." Monatshefte Fur Chemie 125(2): 167-188.
            [2] Agius, P., et al. (2010). "Comparing RNA secondary structures using a relaxed base-pair score." RNA 16(5): 865-878.
        '''
        if method in ['top','rbp']:
            dist = method == 'top' and wRNA.distance or Algorithm.BP_distance
        else:
            dist = Algorithm.BP_distance
        if len(fs.structures)==1:
            return numpy.matrix([0])
        dim = len(fs.structures)
        dismatrix = numpy.array([0]*dim*dim, dtype=numpy.float)
        dismatrix.shape = (dim, dim)
        for i in xrange(dim):
            dismatrix[i,i] = 0
            for j in xrange(i):
                dismatrix[i,j] = dist(fs.structures[i], fs.structures[j])
                dismatrix[j,i] = dismatrix[i,j]
        return dismatrix
    distanceMatrix=staticmethod(distanceMatrix)
    def BP_distance(S1,S2, t = 1):
        '''
        Relaxed base-pair distance between two structures [1].
        Options:
            t: relaxed parameter.
        References:
            [1] Agius, P., et al. (2010). "Comparing RNA secondary structures using a relaxed base-pair score." RNA 16(5): 865-878.
        '''
        # basepairs for S1
        bps1 = []
        m = re.search('[\(\)]+',S1)
        pos= 0 
        while m is not None :            
            bps1.append((pos+m.start()+1,pos+m.end()))
            pos+=m.end()
            m = re.search('[\(\)]+',S1[pos:])
        # basepairs for S2
        bps2 = []
        m = re.search('[\(\)]+',S2)
        pos= 0 
        while m is not None :            
            bps2.append((pos+m.start()+1,pos+m.end()))
            pos+=m.end()
            m = re.search('[\(\)]+',S2[pos:])
        # pairwise distances
        ds = []
        for i,j in bps1:
            ds.append( min([max(abs(i-i1),abs(j-j1)) for i1,j1 in bps2]))
        for i,j in bps2:
            ds.append( min([max(abs(i-i1),abs(j-j1)) for i1,j1 in bps1]))
        # sort ds
        ds = numpy.array(sorted(ds,reverse=True))
        # calcluate min m
        minm = len(ds)
        for m in range(len(ds)+1):
            if numpy.all(ds[m:] <= t*m):
                minm = m
                break
        return minm
    BP_distance=staticmethod(BP_distance)

class Utils(object):
    ''' Utils for RNA structure prediction. '''
    DegreeAlphabet = ['F','K','C', 'R']
    AnyToC = { 'F': lambda x: (x-32.)/1.8, 'K': lambda x: x-273.16, 'R': lambda x: x*1.25 }
    CToAny = { 'F': lambda x: x*1.8+32., 'K': lambda x: x+273.16, 'R': lambda x: x* 0.8 }
    def TempConverter(temp, _from, _to):
        '''
        Convert temperature among different systems.
            C = Celsius
            F = Fahrenheit
            K = Kelvin (absolute temperature)
            R = Reaumur Equals
        '''
        _from = _from.upper()
        _to = _to.upper()
        if _from not in Utils.DegreeAlphabet or _to not in Utils.DegreeAlphabet:
            raise ValueError("Temperature format not supported. Please choose from ['F','C','K', 'R'].")
        if _from == 'C':
            return Utils.CToAny[_to](temp)
        if _to == 'C':
            return Utils.AnyToC[_from](temp)
        t = Utils.AnyToC[_from](temp)
        return Utils.CToAny[_to](t)
    TempConverter=staticmethod(TempConverter)
    def PARSToMfold(constraints):
        ''' Convert PARS constraints to mfold constraints. '''
        mconstraints = []
        for match in re.finditer(r'(\|+)', constraints):
            mconstraints.append('F {0} 0 {1}'.format(match.start()+1, match.end() - match.start()))
        for match in re.finditer(r'(x+)', constraints):
            mconstraints.append('P {0} 0 {1}'.format(match.start()+1, match.end() - match.start()))
        return '\n'.join(mconstraints)
    PARSToMfold=staticmethod(PARSToMfold)
    def PARSToFold(constraints):
        '''
        Convert PARS constraints to Fold (RNAStructure package) constraints.
        Example of Fold contraints:
            DS:
            2
            3
            -1
            SS:
            5
            7
            -1
        '''
        fconstraints = "DS:"
        for match in re.finditer(r'(\|+)', constraints):
            for i in range(match.start(),match.end()):
                fconstraints += "\n{0}".format(i+1)
        fconstraints += "\n-1\nSS:"
        for match in re.finditer(r'(x+)', constraints):
            for i in range(match.start(),match.end()):
                fconstraints += "\n{0}".format(i+1)
        fconstraints +="\n-1\nMod:\n-1\nPairs:\n-1 -1\nFMN:\n-1\nForbids:\n-1 -1"
        return fconstraints
    PARSToFold=staticmethod(PARSToFold)
    def ct2dot(ctname):
        ''' 
        CT format to DOT format.
        Input can be CT file names or a CT file string.
            50      dG = -6.2       seqname
            1       A       0       2       0       1       0       0
            2       U       1       3       0       2       0       0
            3       G       2       4       0       3       0       0
            4       A       3       5       0       4       0       0
            5       C       4       6       0       5       0       0
            6       A       5       7       0       6       0       0
            7       C       6       8       0       7       0       8
            8       A       7       9       36      8       7       9
            9       G       8       10      35      9       8       10
            10      C       9       11      34      10      9       11
            11      U       10      12      33      11      10      12
            12      U       11      13      31      12      11      13
            13      C       12      14      30      13      12      14
            14      A       13      15      29      14      13      15
        Note:
            (1) dG can be ENERGY, energy or dG.
            (2) the first 6 fields are required.
            (3) for RNA without structure predicted, no dG field provided.
        For pseudoknots like this:
            '(((..[[[...)))...]]]'
            The first seen pairs are labelled with '()', and the following conflict ones are labelled as '[]'
        '''
        structures = []
        scores = []
        try: # file name
            fh = open(ctname)
            cts = fh.readlines()
            fh.close()
        except: # string of CTs
            cts = ctname.split('\n')
        A = []
        B = []
        score = 0.
        seq = ''
        name = ''
        # suppose there might be multiple structures in the ct file.
        idx = 0
        while idx < len(cts):
            line = cts[idx].split()
            if len(line)>=6 and line[0] == line[5]:
                if line[5] == '1': # first line
                    if len(A) > 0:
                        s = Utils._ABToDot(A,B)
                        if s.find('(') != -1 or s.find('[') != -1:
                            structures.append(s)
                            scores.append(score)
                    A = [int(line[0])]
                    B = [int(line[4])]
                    seq += line[1]
                    # Parse header                
                    name = cts[idx-1].split()[-1]
                    if name.strip('.') == "":
                        name = "NONAME"
                    m = re.findall('-*\d+\.*\d*',cts[idx-1])
                    if len(m) >= 1:
                        ctlen = int(m[0])
                    if len(m) >=2:
                        score = float(m[1])
                else:
                    A.append(int(line[0]))
                    B.append(int(line[4]))
                    seq += line[1]
            idx +=1
        if len(A) > 0:
            s = Utils._ABToDot(A,B)
            if s.find('(') != -1 or s.find('[') != -1:
                structures.append(s)
                scores.append(score)
        return FastS(name,seq,structures,scores) 
    ct2dot=staticmethod(ct2dot)
    def _ABToDot(A, B):
        ''' Generate DOT structure from a list of pairs. '''
        structure = numpy.repeat('.', len(A))
        bs = [len(A) + 1]
        for a, b in izip(A, B):
            if a > b:
                try:
                    bs.remove(a)
                except:
                    pass
            elif b > min(bs):
                structure[a-1] = '['
                structure[b-1] = ']'
            else:
                structure[a-1] = '('
                structure[b-1] = ')'
                bs.append(b)
        return structure.tostring()
    _ABToDot=staticmethod(_ABToDot)
    def pknots2dot(pkname):
        ''' Convert pknots output to dot format. '''
        structures = []
        scores = []
        try:
            with open(pkname) as fh:
                pks = fh.readlines()
        except:
            pks = pkname.split('\n')
        idx = 0
        while idx < len(pks):
            # Find SEQ
            while idx < len(pks) and not pks[idx].startswith('SEQ'):
                idx += 1
            if idx >= len(pks):
                break
            # Read pairs
            A = []
            B = []
            idx += 1
            while not pks[idx].startswith('---'):
                A.extend([int(i) for i in pks[idx+1].split()])
                B.extend([int(i) for i in pks[idx+2].replace('.','0').split()])
                idx += 4
            # Parse pairs
            structures.append(Utils._ABToDot(A, B))
            # Find energy line
            while not pks[idx].startswith('energy'):
                idx += 1
            scores.append(float(pks[idx].split()[-1]))
        return (structures, scores)
    pknots2dot=staticmethod(pknots2dot)    
    def dot2ct(tfasts): # 
        ''' Convert DOT format to CT format. '''
        ctstring = []
        for st, sc  in izip(tfasts.structures, tfasts.scores):
            # print header
            ctstring.append ("%5d  ENERGY = %-3.2f  %s" % (len(tfasts), sc, tfasts.name))
            stack1=[]
            stack2=[]
            pairs={}
            for i,c in enumerate(st):
                if c == '(':
                    stack1.append(i+1)
                elif c == '[':
                    stack2.append(i+1)
                elif c == ')':
                    pairs[i+1] = stack1.pop()
                    pairs[pairs[i+1]] = i+1
                elif c == ']':
                    pairs[i+1] = stack2.pop()
                    pairs[pairs[i+1]] = i+1
            for i in xrange(1,len(tfasts)+1): # ###24#A######23###25####0###24
                ctstring.append( " %4d %s %7d%4d %4d %4d" % (i, tfasts.seq[i-1], i-1, i+1, pairs.get(i,0),i) )
        return '\n'.join(ctstring)
    dot2ct=staticmethod(dot2ct)
    def draw(fs, ftype = 'ps'):
        '''
        Draw structure in Postscript or SVG format.
        Input: FastS object
            >YIL140W
            AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
            .(((.......)))............(((...)))...............    (-3.1)
            ..............(((.........(((...)))......)))......    (-3.0)
            ...........((((............))))...................    (-3.0)
        Usage:
            Utils.draw(fs, 'png')
        Output:
            YIL140W-1.png  YIL140W-2.png YIL140W-3.png
        '''
        ftype = ftype.lower()
        cts = Utils.dot2ct(fs)
        for i in range(len(fs.structures)):
            fprefix = "{0}-{1}".format(fs.name, i+1)
            if ftype == 'svg':
                wRNA.plot(cts,fprefix+".svg", 1, i+1)
            elif ftype == 'ps':
                wRNA.plot(cts,fprefix+".ps", 0, i+1)
            else:
                lps = wRNA.plot(cts,"", 0, i+1)
                try:
                    p = Popen(['convert', '-', fprefix+"."+ftype], stdin = PIPE, stdout = PIPE, stderr = PIPE)
                    p.stdin.write(lps)
                    pstdout, pstderr = p.communicate()
                except OSError as e:
                    raise OSError(e)
                finally:
                    if pstderr:
                        print >> sys.stderr, "ERROR: convert run error:",pstderr
        return
    draw=staticmethod(draw)

class Test(object):
    ''' Test Module. '''
    def testPredictor():
        ''' Test Predictor. '''
        # test FastC class
        T = 37
        threshold = 0
        tfastc = FastC('YIL140W','AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU','..xx......||..........xx...|......||..............')
        print "RNAfold (Vienna RNA package) Test:"
        print "Sequence for prediction:"
        print tfastc
        
        # test RNAfold
        print "\nRNAfold predicted structure:"
        print Predictor.RNAfold(tfastc, T, threshold)
        print 

        # test Fold
        print "\nFold (RNAStructure package) predicted structure:"
        print Predictor.Fold(tfastc,T, threshold)
        print 

        # test mfold
        print "\nmfold predicted structure:"
        print Predictor.mfold(tfastc,T, threshold)
        print 

        # test UNAfold
        print "\nUNAFold predicted structure:"
        print Predictor.UNAFold(tfastc, T, threshold)
        print 

        # test pknots
        print "\npknots predicted structure:"
        print Predictor.pknots(tfastc, T, threshold)
        print

        # test pknotsRG
        print "\npknotsRG predicted structure:"
        print Predictor.pknotsRG(tfastc, T, threshold)
        print
        
        # test ipknot
        print "\nipknot predicted structure:"
        print Predictor.ipknot(tfastc, T, threshold)
        print 

        # test FoldMerge
        print "\n Merge all the predicted structures together."
        print Predictor.FoldMerge(tfastc, T, predictors = ['RNAfold','pknots','Fold','mfold'])
        print
    testPredictor=staticmethod(testPredictor)
    def testUtils():
        ''' Test Utils. '''
        tfastc = FastC('YIL140W','AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU','.|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........')
        tfasts = Predictor.RNAfold(tfastc, 37, 0)
        print "DOT format to CT format:"
        cts = Utils.dot2ct(tfasts)
        print cts
        print
        print "CT format back to DOT format:"
        print Utils.ct2dot(cts)
    testUtils=staticmethod(testUtils)

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv) == 1:
        sys.exit('''
        Module description:
            RNA Structure Prediction Module.
        Test:
            python wRNA.py test
        
        Usage:
            import wRNA
            tfastc = wRNA.FastC('YIL140W','AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU','.|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........')
            print tfastc
            print wRNA.Predictor.RNAfold(tfastc)

        Output:
            >YIL140W
            AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
            .|xx.|.x..||.|...xx.x.xx...||.x|..||...|..........
            >YIL140W
            AUGACACAGCUUCAGAUUUCAUUAUUGCUGACAGCUACUAUAUCACUACU
            ((...(..((((.((............)))..))))...)).........    (9.70)
        ''')
    else:
        Test.testPredictor()
        #Test.testUtils()

