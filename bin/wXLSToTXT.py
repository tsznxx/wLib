#!/usr/bin/python
#Last-modified: 07 Jan 2014 04:14:37 PM

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

import sys,os
import string,xlrd

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

def basename(filename):
    return os.path.basename(filename).split('.')[0]

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    if len(sys.argv)==1:
        print "Usage: "+sys.argv[0]+" *.xls(x)"
    else:
        book = xlrd.open_workbook(sys.argv[1],formatting_info=False)
        bname=basename(sys.argv[1])
        for sname in book.sheet_names():
            fh=open(bname+"_"+sname.replace(' ','_')+".txt",'w')
            sh = book.sheet_by_name(sname)
            for r in range(sh.nrows):
                fh.write("\t".join([str(elem.value) for elem in sh.row(r)]))
                fh.write("\n")
            fh.close()

