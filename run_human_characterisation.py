#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 11:58:36 2021

@author: rosie
"""
import human_characterisation_functions as hcf


filename='/Volumes/TOSHIBA EXT/Human slices/OP210319/2021_03_19_0014.abf'    
res=hcf.APproperties(filename,4)    

#%%
vctpfile='/Volumes/TOSHIBA EXT/Human slices/OP201027/20o27001.abf'    
channel=4
ResRa, ResRi = hcf.access_resistance(vctpfile, channel)
#%% resting membrane potential

vmfile='/Volumes/TOSHIBA EXT/Human slices/OP201027/20o27001.abf'    
channel=1

rm=hcf.restingmembrane(vmfile,channel)
