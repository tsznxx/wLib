#!/usr/bin/python
#Last-modified: 28 Jan 2014 03:19:42 PM

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

import os,sys
import string
from setuptools import setup, find_packages, Extension

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__ == '__main__':
    if float(sys.version[:3])<2.7 or float(sys.version[:3])>=2.8:
        sys.stderr.write("CRITICAL: Python version must be 2.7!\n")
        sys.exit(1)

    # includepy = "%s/include/python%s" % (sys.prefix, sys.version[:3])
    with open("README",'r') as fh:
        long_description = fh.read()
    # version
    with open("VERSION",'r') as fh:
        version = fh.read().split()[1]

    # Compile Kent lib
    if sys.argv[1] == "install" or sys.argv[1] == 'build':
        print >>sys.stderr, "Compile KentLib ..."
        os.system('cd external/KentLib && make && cd ../..')
        print >>sys.stderr, "Compile Vienna RNA Package ..."
        os.system('cd external/RNAlib && make && cd ../..')
    elif sys.argv[1] == "clean":
        print >>sys.stderr, "Clean KentLib ..."
        os.system('cd external/KentLib && make clean && cd ../..')
        print >>sys.stderr, "Clean Vienna RNA Package ..."
        os.system('cd external/RNAlib && make clean && cd ../..')
        print >>sys.stderr, "Clean dist and egg info ..."
        os.system('if [ -d dist ]; then rm -rf dist; fi')
        os.system('if [ -f wLib.egg-info ]; then rm wLib.egg-info; fi')
        os.system('if [ -d wLib.egg-info ]; then rm -rf wLib.egg-info; fi')
    
    # install requirement
    with open("requirements.txt", "r") as fh:
        install_requires = [x.strip() for x in fh.readlines() if not x.startswith("wLib")]

    setup(name="wLib",
          version=version,
          author='Yunfei Wang',
          author_email='yfwang0405@gmail.com',
          url='http://tsznxx.appspot.com',
          license="UTD",
          keywords = "Python Sequencing Bed",
          description = ("Python Modules for High-throughput Sequencing Data Analysis."),
          long_description = long_description,
          package_dir={'wLib':'src'},
          packages = ['wLib'],
          scripts=['bin/wBamToWig.py',
                   'bin/wBedToWig.py',    
                   'bin/wXLSToTXT.py',
                   'scripts/wBedExtend.py',
                   'scripts/wBedToFasta.py',
                   'bin/wSam2Bed.py',
                   'bin/wWigPlot.py',
                   'scripts/wFindNearestAnnotation.py',
                   'scripts/wGetTSS.py',
                   'scripts/wlncRNA.py',
                   'scripts/wGTFToTab.py',
                   'scripts/wTabToBed.py',
                   'scripts/wFormatFasta.py',
                   'scripts/wGCContent.py',
                   'scripts/wGetSeqByCoordinates.py',
                   'scripts/wGetSeqByName.py',
                   'scripts/wRandomBed.py',
                   'scripts/wFindNearestTwoAnnotation.py',
                   'scripts/wFindNearestAnnotation.py',
                   'scripts/wBedAnnotation.py',
                   'scripts/wBedGetWig.py',
                   'scripts/wGetIntron.py',
                   'scripts/wGetExon.py',
                   'scripts/wMergeSort.py'],
          ext_modules=[Extension('wWigIO',['external/KentLib/wWigIO/wWigIO.c'],extra_link_args=['-DMACHTYPE_x86_64','-lz','-lm','external/KentLib/lib/jkweb.a'],extra_compile_args='-w -shared -fPIC -p -Iexternal/KentLib/inc'.split(' ')),
                       Extension('wTwoBitIO',['external/KentLib/wTwoBitIO/wTwoBitIO.c'],extra_link_args=['-DMACHTYPE_x86_64','-lz','-lm','external/KentLib/lib/jkweb.a'],extra_compile_args='-w -shared -fPIC -p -Iexternal/KentLib/inc'.split(' ')),
                       Extension('wRNA',['external/RNAlib/wRNA/wRNA.cpp'],extra_link_args=['-lm','external/RNAlib/libRNA.a'],extra_compile_args='-w -shared -fPIC -p -Iexternal/RNAlib/fold -Iexternal/RNAlib/plot -Iexternal/RNAlib/wRNA'.split(' '))
                       ],
          classifiers=['Development Status :: 5 - productive',
                       'Environment :: Console',
                       'Intended Audience :: Developers',
                       'License :: OSI Approved :: Artistic License',
                       'Operating System :: CentOS :: CentOS release 6.4 (Final)',
                       'Operating System :: Fedora :: Fedora release 17 (Final)',
                       'Operating System :: Ubuntu :: Ubuntu 12.04',
                       'Operating System :: RedHatEnterpriseServer :: Red Hat Enterprise Linux Server release 5.5 (Tikanga)',
                       'Programming Language :: Python',],
          install_requires=install_requires)

