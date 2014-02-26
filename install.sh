#!/bin/sh
#Last-modified: 26 Feb 2014 11:14:47 AM

####################### Module/Scripts Description ######################
#  
#  Copyright (c) 2008 Yunfei Wang <tszn1984@gmail.com>
#  
#  This code is free software; you can redistribute it and/or modify it
#  under the terms of the BSD License (see the file COPYING included with
#  the distribution).
#  
#  @status:  experimental
#  @version: $Revision$
#  @author:  Yunfei Wang
#  @contact: tszn1984@gmail.com
#
#########################################################################


USAGE=" Usage 1: install to \$HOME/local\n > $0 install\n\nUsage 2: install to a specified path\n > $0 path\n"
case $# in
	0) echo -en $USAGE
	   exit;;
	*) ;;
esac

# install path
install_path=$1
if [ "$install_path" == "install" ]; then
	install_path=$HOME/local
fi

# version
version=$(more VERSION|cut -f 2)
echo $version

# uninstall the old version
cat installed_files.txt|xargs rm -rf

# install current version
python setup.py install --record installed_files.txt --prefix=$install_path
