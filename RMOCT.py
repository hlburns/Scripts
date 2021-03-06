#! /usr/bin/env ipython 
################################################################                                     
##                                                            ##                                     
##      RMOC General purpose script for RMOC calculation      ##                                     
##                        Helen Burns                         ##   
################################################################


''' Please give a run name and the this will then calculate 
    the RMOC in temperature layers and draw a picture for 
    all files in that folder and save the figures in ~/Figures/Folder.                               

    For spinup timeseries please use the Spin up diagnostics notebook.''' 
###################################################################   
###################################################################   
#--Import modules--# 
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob
#--Take terminal inputs--# 
if len(sys.argv) < 2:
   print '''\n    OOPSIE                                                                             
            usage: Go to the folder with netcdf Psi output                                          
            (named *Psi.nc with a grid folder located ../ or   
            within the folder (for working directories) \n'''
   sys.exit(1)
OP = sys.argv[1]
#--Set folder structure--# 
x=os.getcwd()
lists=glob.glob(x+'/'+OP+'/*Psi.nc')
#--Main For loop--#                                                           
#For every .nc file in the folder                                             
#Read in netcdf variables                                                     
#Decide resolution and time step                                              
#Find grid info                                                              
#Calculate the streamfunction                                                 
#Draw a picture and save it
if not lists:
   print '''\n    OOPSIE:                                                                            
            no *Psi.nc folders here!!                                                                
            usage: Go to the folder with netcdf Psi output                                           
            (named *all.nc with a grid folder located ../ or       
            within the folder (for working directories) \n'''
   sys.exit(1)
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]
X=file2.variables['X']
X=X[:]
if len(X)>400 :
    Q_levs = (np.arange(-50,50,2.5))
    Psi_levs = Q_levs / 10
    Q_ticks = np.arange(-50,50,5.)
    Psi_ticks = Q_ticks / 10
else:
    Q_levs = (np.arange(-10,10)+0.5)
    Psi_levs = Q_levs / 10
    Q_ticks = np.arange(-8,9,2.)
    Psi_ticks = Q_ticks / 10
for file in lists:
    file2read = netcdf.NetCDFFile(file,'r')
    lvrho=file2read.variables["LaVH1TH"]
    lvrho=lvrho[:]
    time=file2read.variables['T']
    ti=time[:]
    dx=Y[1]-Y[0] # Find Resolution
    VT=np.sum(lvrho*dx,axis=3) #integrate Vdx along x
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column
    psi=np.mean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv and put back in right order
    y=Y/1000
    z=np.array(range(1,60))
    Rho = np.genfromtxt(x+'/'+OP+'/Temp', delimiter = ',') 
    nolayers=len(psi[:,1])
    Rho=Rho[0:nolayers]#The layers package bins a layer so adjust for that
    start=int(np.divide(ti[0],(86400*360)))#Find run start and stop times
    end=int(np.divide(ti[-1],(86400*360)))
    #if np.max(abs(psi))>3:
    #   Psi_levs = (np.arange(-25,25)+0.5)/10
    cf=plt.contourf(y,Rho,psi,Psi_levs,cmap=plt.cm.seismic) #Use b2r colourmap 
    #plt.clim(-1,1) # Put 0 to white
    cbar=plt.colorbar(cf, ticks=Psi_ticks, shrink=0.8)
    plt.title("RMOC for years "+str(start)+"-"+str(end)+" ")
    plt.xlabel('Distance (km)')
    plt.ylabel('Temperature $^o$C')
    cbar.ax.set_ylabel('$\psi \,\, (sv)$')
    figpath=( os.path.expanduser('~')+"/Figures/"+OP)
    if not os.path.exists(figpath):
          os.makedirs(figpath)
    plt.savefig(figpath+"/Psires"+str(start)+"-"+str(end)+".png")
    plt.clf()
    

