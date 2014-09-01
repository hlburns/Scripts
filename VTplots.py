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
#import csv
import sys
from pylab import *
from numba import autojit
#import glob
###
OP="3daynokpp"
OP2 = sys.argv[1]
Comp='Mobilis'
x=os.getcwd()
#lists=glob.glob(x+'/'+str(OP)+'/*all.nc')
###
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Z=file2.variables['Z']
Yc=file2.variables['Y']
Yc=Yc[:]*1
y=Yc/1000
z=Z[:]*1
#Read in the VTtav files
fileVTp = netcdf.NetCDFFile(x+'/'+str(OP)+"/VTprimebar.nc",'r') 
VTbar=fileVTp.variables['VT']
VTzone=np.mean(VTbar[:]*1,axis=2)
cf1=plt.contourf(y,z,VTzone,100,cmap=plt.cm.seismic)
clim(-0.1,0.1) # Put 0 to white                                                                   
cbar = plt.colorbar(cf1, shrink=0.8)
plt.title("V'T' ("+OP+")")
plt.ylabel("Depth")
plt.xlabel("Meridional Distance (km)")
x=( os.path.expanduser('~')+"/Figures/"+Comp+"/"+OP)
if not os.path.exists(x):
          os.makedirs(x)
plt.savefig(x+'/'+Comp+"/"+OP+"/VTplt.png")
#Read in the VTtav files
fileVT = netcdf.NetCDFFile(x+'/'+str(OP2)+"/VTprimebar.nc",'r') 
VTbar=fileVT.variables['VT']
VTzone2=np.mean(VTbar[:]*1,axis=2)
Anom=abs(VTzone2)-abs(VTzone)
Z=Z[:]*1
cf1=plt.contourf(y,Z,Anom,40,cmap=plt.cm.seismic)
clim(-0.025,0.025) # Put 0 to white                                                                  
cbar = plt.colorbar(cf1, shrink=0.8)
plt.title("V'T' anomally ("+OP2+")")
plt.ylabel("Depth")
plt.xlabel("Meridional Distance (km)")
x=( os.path.expanduser('~')+"/Figures/"+Comp+"/"+OP2)
plt.savefig(x+"/VTanomplt.png")
