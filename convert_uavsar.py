#!/Library/Frameworks/Python.framework/Versions/Current/bin/python

###############################################################################
# convert_uavsar.py 
#
# Project:   
# Purpose:  Converts UAVSAR .grd files into usable GeoTiffs with speckle filtering 
# Author:   Scott Arko
#
###############################################################################
# Copyright (c) 2010, Scott Arko 
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Library General Public License for more details.
# 
# You should have received a copy of the GNU Library General Public
# License along with this library; if not, write to the
# Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# Boston, MA 02111-1307, USA.
###############################################################################
# Notes:
#
# This program makes extensive use of ideas, algorithms and sometimes direct code
# copy from gdal_merge.py as distributed with the gdal suite
# 
# Versions:
# 
# Dec 2010 -- Original creation (v0.1), S. Arko
#       This version reads the whole grd file at one time.  A potential problem exists where files are
#       just too big to be loaded into memory.  Next version will incorporate reading by line instead of 
#       whole images


#####################
#
# Import all needed modules right away
#
#####################

from osgeo import gdal
import numpy as np
import os
import sys
import math
import saa_func_lib as sa


def Usage():
        print '*******************************'
        print '*** USAGE ******'
        print '*** convert_uavsar.py ann_file grd_file <-filter>'
        print '-filter option applies a two-dimmensional modified boxcar filter to the data prior to exporting'
        print '*******************************'

if len(sys.argv) <3:
        Usage()
        sys.exit( 0 )

# Get the annotation file and data file from the command line
annfile = sys.argv[1]
file = sys.argv[2]

argv = gdal.GeneralCmdLineProcessor( sys.argv )

# Parse command line arguments.
i = 3
filter = 0
tif = 0
while i < len(argv):
    arg = argv[i]
    if arg == '-o':
        i = i + 1
        out_file = argv[i]
        
    # The -ot option is not implemented at this time, but would be used for specifying other output 
    # formats if desired
    elif arg == '-ot':
        i = i + 1
        band_type = gdal.GetDataTypeByName( argv[i] )
        if band_type == gdal.GDT_Unknown:
                print 'Unknown GDAL data type: ', argv[i]
                sys.exit( 1 )

    elif arg == '-filter':
        i = i+1
        filter = 1

    elif arg == '-tif':
        i = i + 1
        tif = 1 

    elif arg[:1] == '-':
        print 'Unrecognised command option: ', arg
        Usage()
        sys.exit( 1 )

# parse text file

print 'Parsing text file'
logfile = open(annfile, 'r').readlines()

if tif==0:
    grdfile = file
else:
    grdfile = file.replace('tif','grd')

# Gather information from the annotation file
KEYWORDS = ['grd_pwr.set_rows', 'grd_pwr.set_cols','grd_pwr.row_addr','grd_pwr.col_addr','grd_pwr.row_mult','grd_pwr.col_mult',grdfile]
for line in logfile:
    temp = line.split()
    if len(temp) > 0:
        if temp[0] in KEYWORDS[0]:
                rows = temp[3]
        if temp[0] in KEYWORDS[1]:
                cols = temp[3]
        if temp[0] in KEYWORDS[2]:
                ul_lat = float(temp[3])
        if temp[0] in KEYWORDS[3]:
                ul_lon = float(temp[3])
        if temp[0] in KEYWORDS[4]:
                lat_step = float(temp[3])
        if temp[0] in KEYWORDS[5]:
                lon_step = float(temp[3])
    if len(temp) > 4:
        if temp[2] == KEYWORDS[6]:
            data_size = temp[6]

print 'Data size is ',data_size
# Figure out data type data can either be float*32 or double*64

word_size = (float(data_size)/(float(rows)*float(cols)))

if(word_size % 2 == 0):
    word_size = int(word_size)
    if word_size == 4:
        data_type = 4
        bands = 1
        band_type = 'bsq'
    elif word_size == 8:
        data_type = 4
        bands = 2
        band_type = 'bip'
else:
    print 'Cannot resolve the word size to an integer'
    sys.exit(1)

# Create ENVI-style header file to allow GDAL to read the raw data
# as some format that is known to it.  This is a bit of a cheat, but 
# is a really easy way to read straight binary data that you don't know a 
# whole lot about

if (tif == 0):
        print 'Creating header file'
        hdrfile = file.replace('grd','hdr')
        hdr = open(hdrfile,'w')
        hdr.write('ENVI\n')
        hdr.write('description = {'+file+'}\n')
        hdr.write('samples = '+str(cols)+'\n')
        hdr.write('lines = '+str(rows)+'\n')
        hdr.write('bands   = '+str(bands)+'\n')
        hdr.write('header offset = 0\n')
        hdr.write('file type = ENVI Standard\n')
        hdr.write('data type = '+str(data_type)+'\n')
        hdr.write('interleave = '+band_type+'\n')
        hdr.write('byte order = 0\n')
        hdr.write('band names = {\n')
        hdr.write('}\n')
        hdr.write('\n')
        hdr.close()
        
# Open grd file and read data
print 'Processing files ',file
oh = sa.open_gdal_file(file)

if bands == 1:
        (ox,oy,proj,trans,odata) = sa.read_gdal_file(oh,1)
        np.putmask(odata,odata>5,0)
elif bands ==2:
        (ox,oy,proj,trans,odata) = sa.read_gdal_file(oh,1)
        np.putmask(odata,odata>5,0)
        (ox,oy,proj,trans,pdata) = sa.read_gdal_file(oh,2)

# Create geotransform based on ann file and simple EPSG:4326 WKT for projection
otrans = [ul_lon,lon_step,0,ul_lat,0,lat_step]
oproj = 'GEOGCS[\"WGS 84\",DATUM[\"WGS_1984\",SPHEROID[\"WGS 84\",6378137,298.257223563,AUTHORITY[\"EPSG\",\"7030\"]],AUTHORITY[\"EPSG\",\"6326\"]],PRIMEM[\"Greenwich\",0,AUTHORITY[\"EPSG\",\"8901\"]],UNIT[\"degree\",0.01745329251994328,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4326\"]]'

outfile = file.replace('grd','tif')

if (tif == 1):
    outfile = file.replace('tif','filt.tif')

if filter==1:
        print 'Applying boxcar filter in x direction'
        outdata = sa.boxcar_y(odata,3)
        print 'Applying boxcar filter in y direction'
        outdata2 = sa.boxcar_x(odata,3)
        out = (outdata + outdata2)/2
        sa.write_gdal_file_float(outfile,otrans,oproj,out)
else:
        sa.write_gdal_file_float(outfile,otrans,oproj,odata)

if bands == 2:
        outfile = file.replace('grd','phase.tif')
        sa.write_gdal_file_float(outfile,otrans,oproj,pdata)
