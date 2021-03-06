#! /usr/bin/env ipython
################################################################
##                                                            ##
##     Create time sereis of Av EKE                           ##
##                      Helen Burns                           ##
################################################################


''' Please give a run name in order to file away the output in
     ~/Figures/Runname
    This script will then calculate some spinup diagnostics 
    and draw a picture for all files in
    that folder and save the figures in ~/Figures/Folder.'''    
###################################################################        
################################################################### 
#--Import modules--#
from scipy.io import netcdf
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
import matplotlib.pyplot as plt
import os
from numba import autojit
import csv
import sys
import glob
#--Take terminal inputs--#
if len(sys.argv) < 2:
   print '''\n    OOPSIE
            usage: Go to the folder with netcdf V output
            (named *all.nc with a grid folder located ../ or
            within the folder (for working directories) \n'''
   sys.exit(1)
OP = sys.argv[1]
#--Set folder structure--#
x=os.getcwd()
lists=glob.glob(x+'/*all.nc')
#--Main For loop--#
#For every .nc file in the folder
#Read in netcdf variables
#Decide resolution
#Find grid info
#Calculate the streamfunction
#Draw a picture and save it
if not lists:
   print '''\n    OOPSIE:
            no *all.nc folders here!!
            usage: Go to the folder with netcdf V output
            (named *all.nc with a grid folder located 
            within the folder (for working directories) \n'''
   sys.exit(1)
file2=netcdf.netcdf_file('grid.nc','r')
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
file2=netcdf.netcdf_file('grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]
def regridy(Variable):
    Vc=(Variable[:,:,0:-1,:]+Variable[:,:,1::,:])/2
    return Vc
def regridx(Variable):
    Vc=(Variable[:,:,:,0:-1]+Variable[:,:,:,1::])/2
    return Vc
numba_regridy = autojit()(regridy)
numba_regridy.func_name = "numba_regridy"
numba_regridx = autojit()(regridx)
numba_regridx.func_name = "numba_regridx"
## EKE ##
def EKEf(file):
       file2read = netcdf.NetCDFFile(file,'r')
       V=file2read.variables['VVEL']
       U=file2read.variables['UVEL']
       V=V[:]*1
       U=U[:]*1
       Vc=numba_regridy(V)
       Uc=numba_regridx(U)
       #Find timestep and resolution
       time=file2read.variables['T']
       ti=time[:]
       dx=Y[1]-Y[0]
       k=0.5*(Uc**2+Vc**2)
       EKE=np.nanmean(k)
       T=ti[0]
       start=int(np.divide(ti[0],(86400*360)))
       end=int(np.divide(ti[-1],(86400*360)))
       return EKE,T

EKEts=[]
TS=[]
for file in lists:
       EKE,T=EKEf(file)
       EKEts.append(EKE)
       TS.append(T)
EKEts=np.array(EKEts)
TS=np.array(TS)
year=TS/(86400*360)
a=np.array((year,EKEts))
b=a[:,np.argsort(a[0,:])]
plt.plot(b[0,:],b[1,:])
plt.title("Av domain EKE for years "+str(int(b[0,0]))+"-"+str(int(b[0,-1])))
plt.xlabel('Years')
plt.ylabel('EKE')
x=( os.path.expanduser('~')+"/Figures/"+OP)
if not os.path.exists(x):
          os.makedirs(x)
plt.savefig(x+'/EKEtimeseries.png')
plt.clf()
