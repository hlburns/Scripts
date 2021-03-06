#! /usr/bin/env ipython 
################################################################                                     
##                                                            ##                                     
##      Eddy General purpose script for MOC calculation       ##                                     
##                        Helen Burns                         ##   
################################################################
## !!!!PLEASE NOTE!!!! ##
## YOU MUST USE SciPy 14.0 to run this script ##

''' Please give a run name and the this will then calculate 
    the Eddy MOC and draw a picture for 
    all files in that folder and save the figures in ~/Figures/Folder.                               
    For spinup timeseries please use the Spin up diagnostics notebook.''' 
###################################################################   
###################################################################   
#--Import modules--# 
from scipy.io import netcdf
import numpy as np
from scipy.interpolate import interp1d
from scipy import interpolate
from numba import autojit
import matplotlib.pyplot as plt
import os
import csv
import sys
import glob
from pylab import *
#--Take terminal inputs--# 
if len(sys.argv) < 2:
   print '''\n    OOPSIE                                                                             
            usage: *all.nc *Psi.nc grid.nc Temp must be here'''
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
            Wrong folder you ninny!*all.nc *Psi.nc grid.nc Temp must be here'''
   sys.exit(1)
def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx] 
def regrid(Variable):
    Vc=(Variable[:,0:-1]+Variable[:,1::])/2
    return Vc
numba_regrid = autojit()(regrid)
numba_regrid.func_name = "numba_regrid"
file2=netcdf.netcdf_file(x+'/'+str(OP)+'/grid.nc','r')
Z=file2.variables['Z']
Z=Z[:]*1
Zp1=file2.variables['Zp1']
Zp=Zp1[:]
Y=file2.variables['Yp1']
Y=Y[:]*1
X=file2.variables['X']
X=X[:]*1
lm=file2.variables['HFacC']
Yc=file2.variables['Y']
Yc=Yc[:]
Z=file2.variables['Z']
Z=Z[:]*1
X=file2.variables['X']
X=X[:]*1
lm=lm[:]
lmc=np.array(lm)
lmc[lmc<1]=np.nan
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
    #Stop memory issues!!
    if len(ti)>300:
       continue
    dx=Y[1]-Y[0] # Find Resolution
    VT=np.sum(lvrho*dx,axis=3) #integrate Vdx along x
    VTfdz=np.cumsum(VT[:,::-1,:],axis=1) #sum up the water column
    psi=np.mean(VTfdz[:,::-1,:],axis=0)/10**6 #Time average and put into Sv and put back in right order
    y=Y/1000
    Rho = np.genfromtxt(x+'/'+str(OP)+'/Temp', delimiter = ',') 
    nolayers=len(psi[:,1])
    Rho=Rho[0:nolayers]#The layers package bins a layer so adjust for that
    start=int(np.divide(ti[0],(86400*360)))#Find run start and stop times
    end=int(np.divide(ti[-1],(86400*360)))
    #Now find the corresponding all file! (turn to list and replace Psi with all!)
    a=list(file)
    a[-6]='a'
    a[-5]='l'
    a[-4]='l'
    A=''.join(a)
    if not os.path.isfile(A):
       print "Warning: file "+A+" not found, so i'm skipping this!"
       continue
    file2 = netcdf.NetCDFFile(A,'r')
    Temp=file2.variables['THETA']
    Temp=Temp[:]
    Tav=np.mean(Temp,axis=0)
    V=file2.variables['VVEL']
    V=V[:]*1
    Tavlat=np.mean(Tav,axis=2)
    Yc=Yc/1000
    lmav=np.mean(lmc,axis=2)
    Vtave=np.mean(V,axis = 0)
    #Vtave[Vtave==0]=np.nan
    Vzone=np.nansum(Vtave*dx,axis = 2)
    dz=Zp[0:len(Zp)-1]-Zp[1:len(Zp)]
    # No more super slow forloop!
    psi2=np.apply_along_axis(np.multiply,0,Vzone,dz)
    psi3=np.cumsum(-psi2[::-1,:],axis=0)
    npad = ((0,1), (0,0))
    psi4 = np.pad(psi3, pad_width=npad, mode='constant', constant_values=0)
    y =Y/1000
    Psi=psi  #Convert to Sv 
    Psi2=psi4/10**6
    Psi=numba_regrid(Psi)
    Psi2=numba_regrid(Psi2)
    #Remap Part
    #Expand temperature co-ordinates (30 lvls to 168 lvls)
    Z2=interp1d(Z,Z,axis=0)
    Znew=linspace(int(Z[0]),int(Z[-1]),168)
    Zexp=Z2(Znew)
    T2=interp1d(Z,Tavlat,axis=0)
    Tnew=linspace(int(Z[0]),int(Z[-1]),168)
    Texp=T2(Tnew)
    R2=interp1d(Rho,Rho,axis=0)
    Rnew=linspace(Rho[0],Rho[-1],168)
    Rexp=R2(Rnew)
    P2=interp1d(Rho,psi,axis=0)
    Pnew=linspace(Rho[0],Rho[-1],168)
    Pexp=P2(Pnew)
    Psimap=np.zeros(shape(Texp))
    for i in range(len(Yc)):
        for k in range(len(Zexp)):
            D=Texp[k,i]
            if np.isnan(D):
               Psimap[k,i]=np.nan
            else:
                P=Pexp[:,i]
                I=find_nearest(Rexp, D)
                b=nonzero(Rexp==I)[0][0]
                Psimap[k,i]=P[b]
    #Psimapped=np.multiply(Psimap,lmav) HFacS is req for topo runs
    Psimap[Psimap==0]=np.nan
    #Now put the MOC into 168 lvls to make the eddy plot
    P3=interp1d(Zp,Psi2,axis=0)
    P3new=linspace(Z[0],Zp[-1],168)
    P3exp=P3(P3new)
    Psied=Psimap-P3exp
    x2=( os.path.expanduser('~')+"/Figures/"+OP)
    if not os.path.exists(x2):
          os.makedirs(x2)
    if np.max(abs(psi))>1:
       Psi_levs = (np.arange(-25,25)+0.5)/10
    cf1=plt.contourf(Yc,Zexp,Psimap,Psi_levs,cmap=cm.seismic) #Use b2r colourmap                     
    #clim(-2.5,2.5) # Put 0 to white
    cbar = plt.colorbar(cf1, ticks=Psi_ticks, shrink=0.8)                           
    title("ROC (y,z) year "+str(start)+" "+OP)
    xlabel('Distance (km)')
    ylabel('Depth (m)')
    cbar.ax.set_ylabel('Psi (sv)')
    plt.savefig(x2+"/ROCremap"+str(start)+"-"+str(end)+".png")
    plt.clf()
    cf=plt.contourf(Yc,Zexp,Psimap,Psi_levs,cmap=cm.seismic)
    cbar = plt.colorbar(cf, ticks=Psi_ticks, shrink=0.8)
    clim(-1,1)
    plt.title("Remapped Stream Function "+str(start)+"yrs  "+OP)
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    cbar.ax.set_ylabel('$\psi_{res}$ (Sv)')
    plt.savefig(x2+"/Psiremapped"+str(start)+"-"+str(end)+".png")
    plt.clf()    

