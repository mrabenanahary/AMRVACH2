#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 20:04:46 2022

@author: mrabenanahary
"""

import io
import sys
import os
import numpy as np

# sys.argv[1] must be $AMRVAC_DIR
# sys.argv[2] must be the name of the .par
AMRVAC_DIR = str(sys.argv[1])
ROOT_DIR = '/mysimus/HD/HH_YSO_Jet'
output_dir = 'Output/'
PAR_DIR = ROOT_DIR+'/'+str(sys.argv[2])+str(sys.argv[3])+'.par'
outnum_density = int(sys.argv[4])
outnum_temperature = int(sys.argv[5])
namedir = str(sys.argv[6])
iline_of_output = 15
iline_of_restart = 19
iline_of_timemax = 59
iline_of_timedt = 46
iline_of_density = 303
iline_of_temperature = 313


temperature = np.logspace(1, 9, 16)

density_zero = 2.2008179078168482e-24

nH_log_factor = np.arange(0,9)

density_list = density_zero*10**nH_log_factor

#print(AMRVAC_DIR)


# with is like your try .. finally block in this case
with open(AMRVAC_DIR+PAR_DIR, 'r') as file:
    # read a list of lines into data
    data = file.readlines()

#[print(str(i)+' '+el) for i,el in enumerate(data)]

# now change the 2nd line, note that you have to add a newline
#data[1] = 'Mage\n'

#[print(data[el].split()) for i,el in enumerate(ilines_to_modify)]

#timemax = float((((data[iline_of_timemax].split('='))[1]).split('!')[0]).replace('d','e'))
#print(timemax)

#timedt = float((((data[iline_of_timedt].split('='))[1]).split('!')[0]).replace('d','e'))
#print(timedt)



#data[iline_of_timemax] = 'time_max = '+ str(float((outnum+1)*timedt))+'\n'

filename_suffix = ((data[iline_of_output].split('=')[1]).split('/')[-1]).split("'")[0]
output_dir


filename_suffix_to_restart = filename_suffix+"_nH_{:04d}".format(nH_log_factor[outnum_density])
filename_suffix_to_restart += "_T_{:04d}".format(outnum_temperature)+"_"

filename = namedir+'.txt'

data[iline_of_output] = 'base_filename = \''+output_dir+"nH_{:g}".format(nH_log_factor[outnum_density])+'/'+filename_suffix+"T_{:04d}".format(outnum_temperature)+"_"'\'\n'

#filename_suffix_to_restart = filename_suffix+"{:04d}".format(outnum)+".dat"
#restart_strg = data[iline_of_restart].split('!')[1]
#print(restart_strg,filename_suffix_to_restart)

density_strg = [("{:.8e}".format(density_list[i])).replace('e','d') for i in range(np.shape(density_list)[0])]

temperature_strg = ["{:.8f}d0".format(temperature[i]) for i in range(np.shape(temperature)[0])]

data[iline_of_density] = 'ism_config%density = '+ density_strg[outnum_density]+'\n'

data[iline_of_temperature] = 'ism_config%temperature = '+ temperature_strg[outnum_temperature]+'\n'

print("Launching run for rho = {:s} g.cm-3, T = {:s} K. Filename : {:s}".format(density_strg[outnum_density],
                                                                                temperature_strg[outnum_temperature],
                                                                                filename))

#if(outnum>0):
#    data[iline_of_restart]=restart_strg[:-2]+filename_suffix_to_restart+"'\n"

# and write everything back
with open(filename, 'w') as file:
    file.writelines(data)
