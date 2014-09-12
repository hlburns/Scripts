#! /usr/bin/env ipython                                                                              
################################################################                                     
##                                                            ##                                     
##        Script to create V'T' time average variable         ##                                    
##                        Helen Burns                         ##                                     
################################################################                                    
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import sys
from pylab import *
from numba import autojit
import glob
###
OP = sys.argv[1]
x=os.getcwd()
lists=glob.glob(x+'/'+str(OP)+'/*all.nc')
###
#Read in the timeaveraged files
file2 = netcdf.NetCDFFile(x+'/'+str(OP)+"/Tav.nc",'r') 
Temp=file2.variables['THETA']
Yc=file2.variables['Y']
Z=file2.variables['diag_levels']
X=file2.variables['X']
V=file2.variables['VVEL']
Yp1=file2.variables['Yp1']
#Regrid V to center points
def regrid(Variable):
    Vc=(Variable[:,:,0:-1,:]+Variable[:,:,1::,:])/2
    return Vc
numba_regrid = autojit()(regrid)
numba_regrid.func_name = "numba_regrid"
VC=numba_regrid(V[:]*1)
# Remove excess time dimension
VTAVC=np.squeeze(VC)
TTAV=np.squeeze(Temp[:]*1)
# Load in each file and find V'T' and timeaverage it!
lists=glob.glob(x+'/'+str(OP)+'/*all.nc')
VTprimetav20=0
VTbar20=0
if '240-260all.nc' in lists:
    total=len(lists)-1
else:
    total=len(lists)
for file in lists:
    if file=='240-260all.nc':
        continue
    file2 = netcdf.NetCDFFile(file,'r') 
    Temp=file2.variables['THETA']
    V=file2.variables['VVEL']
    Vc=numba_regrid(V)
    Vprime=Vc-VTAVC
    Tprime=Temp[:]*1-TTAV[:]*1
    VTprime=Vprime*Tprime
    VTprimetav=np.mean(VTprime,axis=0)
    VTprimetav20=VTprimetav20+VTprimetav/total
    VTbar=np.mean(Vc*Temp[:]*1,axis=0)
    VTbar20=VTbar20+VTbar/total
# Write to nc format
f=netcdf.netcdf_file(x+'/'+str(OP)+'/VTprimebar.nc','w')
f.createDimension('X',len(VTprimetav20[1,1,:]))
f.createDimension('Y',len(VTprimetav20[1,:,1]))
f.createDimension('Z',len(VTprimetav20[:,1,1]))
VT=f.createVariable('VT','double',('Z','Y','X'))
VT[:]=VTprimetav20
f.close()
# Write to nc format                                                                                  
f=netcdf.netcdf_file(x+'/'+str(OP)+'/VTbar.nc','w')
f.createDimension('X',len(VTbar20[1,1,:]))
f.createDimension('Y',len(VTbar[1,:,1]))
f.createDimension('Z',len(VTbar[:,1,1]))
VT=f.createVariable('VT','double',('Z','Y','X'))
VT[:]=VTbar20
f.close()
