#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 20:04:46 2022

@author: mrabenanahary
"""

import io
import sys
import os

# sys.argv[1] must be $AMRVAC_DIR
# sys.argv[2] must be the name of the .par
AMRVAC_DIR = str(sys.argv[1])
ROOT_DIR = '/mysimus/HD/HH_YSO_Jet'
PAR_DIR = ROOT_DIR+'/'+str(sys.argv[2])+str(sys.argv[3])+'.par'
outnum = int(sys.argv[4])
iline_of_output = 15
iline_of_restart = 19
iline_of_timemax = 59
iline_of_timedt = 46
ilines_to_modify = [iline_of_output,iline_of_restart,iline_of_timemax,iline_of_timedt] 

#print(AMRVAC_DIR)


# with is like your try .. finally block in this case
with open(AMRVAC_DIR+PAR_DIR, 'r') as file:
    # read a list of lines into data
    data = file.readlines()

#[print(str(i)+' '+el) for i,el in enumerate(data)]

# now change the 2nd line, note that you have to add a newline
#data[1] = 'Mage\n'

#[print(data[el].split()) for i,el in enumerate(ilines_to_modify)]

timemax = float((((data[iline_of_timemax].split('='))[1]).split('!')[0]).replace('d','e'))
#print(timemax)

timedt = float((((data[iline_of_timedt].split('='))[1]).split('!')[0]).replace('d','e'))
#print(timedt)

data[iline_of_timemax] = 'time_max = '+ str(float((outnum+1)*timedt))+'\n'

filename_suffix = ((data[iline_of_output].split('=')[1]).split('/')[-1]).split("'")[0]
filename_suffix_to_restart = filename_suffix+"{:04d}".format(outnum)+".dat"
restart_strg = data[iline_of_restart].split('!')[1]
#print(restart_strg,filename_suffix_to_restart)

if(outnum>0):
    data[iline_of_restart]=restart_strg[:-2]+filename_suffix_to_restart+"'\n"

# and write everything back
with open('split_'+str(outnum)+'.txt', 'w') as file:
    file.writelines( data )
