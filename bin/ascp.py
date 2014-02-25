#!/usr/bin/python
#Last-modified: 06 Jan 2012 04:03:31 PM

""" Module/Scripts Description

Copyright (c) 2008 Yunfei Wang <tszn1984@gmail.com>

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file COPYING included with
the distribution).

@status:  experimental
@version: $Revision$
@author:  Yunfei Wang
@contact: tszn1984@gmail.com
"""

# ------------------------------------
# python modules
# ------------------------------------

import os,sys
import string

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

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" file1 ./"
        print "file:  anonftp@ftp-trace.ncbi.nlm.nih.gov:file1 "
    else:
        srapath=sys.argv[1]
        try:
            place=sys.argv[2]
        except:
            place='./'
        srapath=srapath.split('/')[1:]
        srapath='/'.join(srapath)
        srapath="~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.putty -k 1 -QTr -l200m "+"anonftp@ftp-trace.ncbi.nlm.nih.gov:"+srapath+" "+place
        print srapath
        os.system(srapath)

