#! /usr/bin/env ipython 
################################################################
##                                                            ##
##       Make a movie of temperature evolution with time!     ##      
##                      Helen Burns                           ##
################################################################


''' Give this script the name of your run and run it in the parent dir
    It will go through and open all the *all.nc files and draw pictures
    of each snap shot of the temp field. I've set this to make a movie 
    folder ~/Movies/OP and put the jpg there numbered 1-3000 or whatever!
'''                                           
###################################################################                                  
###################################################################                                  
#--Import modules--#
from scipy.io import netcdf
import numpy as np
from pylab import clf, plot, show, floor, ceil, imshow
from matplotlib.pyplot import xlabel, ylabel, legend, savefig, colorbar, title, clim, pcolor, cm, contourf
import os
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
lists=glob.glob(x+'/'+str(OP)+'/*all.nc')
moviefolder=( os.path.expanduser('~')+"/Movies/Temp"+OP)
if not os.path.exists(moviefolder):
    os.makedirs(moviefolder)
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
            (named *all.nc with a grid folder located ../ or                                                                                
            within the folder (for working directories) \n'''
   sys.exit(1)
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]*1
Y=file2.variables['Y']
Y=Y[:]
for file in lists:
       file2read = netcdf.NetCDFFile(file,'r')
       Temp=file2read.variables['THETA']
       Temp=Temp[:]
       #Find timestep and resolution
       time=file2read.variables['T']
       ti=time[:]
       Tzone=np.nanmean(Temp,axis = 3)
       # Start off with a forloop but there must be a better way!!    
       y =Y/1000
       Q_levs = np.arange(0,8,0.25)
       Q_ticks = np.arange(0,8,1)
       for i in range(len(ti)-1):# bin the last one it's going to be repeated!
           num=ti[i]/(86400*30)#Month number
           cf=contourf(y,Z,Tzone[i,:,:],Q_levs,cmap=cm.seismic) #Use b2r colourmap
           cbar = colorbar(cf, ticks=Q_ticks, shrink=0.8)
           title("Temp field month "+str(num))
           xlabel('Distance (km)')
           ylabel('Depth (m)')
           cbar.ax.set_ylabel('$Temp ^oC$')
           savefig(moviefolder+'/'+str(num)+'.jpg')


           clf()
       
