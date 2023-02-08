#!/bin/python
#
# Product: EARS (Earthscope Automated Receiver Survey)
# 
# Script: getEarsRFByStation.py
#
# Description: script for batch download of the receiver frunctions from EARS for the given station
#
# Author: H. Philip Crotwell, University of South Carolina
#
# Usage: 1. on the EARS main page, search for the desired stations using the search form on the 
#          left panel. and save the "network,station" pairs as "stationList.csv" file.
#
#        2. run this python script in the same directory as the stationList.csv file. 
#           The script will print one wget command per event to download a zip file of 
#           all receiver functions . 
#
#        3. create a shell script of the wget commands and run
#
# Questions: product@iris.washington.edu
#
# History:
#      2011-01-24: released


#!/usr/bin/env python
#%%
import os
import csv
import shutil
import time 
# from comebineRFTraces import combineRFTraces
import comebineRFTraces 
import warnings

skipIfPresent = True
tracePath = 'Ears/gauss_2.5/' # Path where we should find net.sta/traces after extracting. e.g. Ears/gauss_2.5/US.CEH/trace_stuff


"""
 Name: getEarsRFByStation.py - is a Python script  for batch download of the receiver functions from EARS for the
 given station.

 Copyright (C) 2019  Product Team, IRIS Data Management Center

    This is a free software; you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation; either version 3 of the
    License, or (at your option) any later version.

    This script is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License (GNU-LGPL) for more details.  The
    GNU-LGPL and further information can be found here:
    http://www.gnu.org/

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 Usage: 1. on the EARS main page, http://ears.iris.washington.edu/,
           search for the desired stations using the search form on the 
           left panel. and save the "network,station" pairs as "stationList.csv" file.

        2. run this python script in the same directory as the stationList.csv file.
           The script will print one wget command per event to download a zip file of
           all receiver functions .

        3. create a shell script of the wget commands and run

 Questions: product@iris.washington.edu

 History:
    2022-02-25: Modified by Brennan to download traces 
      (instead of just output the commands to download them), extract zip files, and not
      redownload traces that have already been downloaded. 
    2019-09-12 IRIS DMC Product Team (Manoch): Python 3 release
    2011-01-24 Philip Crotwell, University of South Carolina: initial release

"""

results = csv.reader(open('stationList.csv', 'r'))
data_dir = os.path.join('data', 'RF')

mainDir = os.getcwd()

failSta = []

for row in results:

   if row[0].strip().startswith('#'):
      continue

   if len(row) == 1: # Wrong deliminter. Try splitting into [net, sta]
      row = row[0].split()

   (net, sta) = (row[0], row[1])
   out_dir = os.path.join(data_dir, f'{net}_{sta}')
   out_file = os.path.join(out_dir, f'{net}_{sta}.zip')
   tracePathNetSta = tracePath + '{}.{}'.format(net, sta)

   # Build commands. 
   make_dir = (f'mkdir -p {out_dir}')
   download_data = (f'wget -O {out_file} "http://ears.iris.washington.edu/receiverFunction.zip?netCode={net}&stacode={sta}'
            f'&minPercentMatch=80&sgaussian=2.5"')

   #Check if data is already downloaded. Possibly skip if it is downloaded. 
   if skipIfPresent and os.path.exists(tracePathNetSta + '/rfArr.mat' ):
       print('Already downloaded receiver functions for {}.{}'.format(net, sta))
       continue
   
   #Execute the download commands. 

   # os.popen(make_dir).read()
   # os.popen(download_data).read()
   # shutil.unpack_archive(out_file)
   # os.remove(out_file)

   # tries = 0
   # succeeded = False
   # while not succeeded: 
   #    if tries > 10: 
   #       continue
   #    else:

   try: 
      matArray = comebineRFTraces.combineRFTraces(trPath = tracePathNetSta + '/')
      print('Finished making Matlab structure for {}.{}'.format(net, sta))
   except: 
      warnings.warn("""brb2022.02.27 There was a glitch (probably Segmenation fault 11?) in combineRFTraces. 
                  You Probably don't have matlab receiver function file. Sometimes trying again works.""")
      failSta.append(row)

if len(failSta) > 0: 
   print('Failures: ') 
   print(failSta) 


# %%
