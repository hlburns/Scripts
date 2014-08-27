#! /usr/bin/env ipython  
''' This Script is calculate the maximum overturning streamfunctions for each model set up'''
from scipy.io import netcdf
from numba import autojit
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys
##~~~~~~~~~~Set up and sanity checks~~~~~~~~~~~##
sys.path.append(os.path.expanduser('~')+"/bin/")
from Calculations import MOC, find_nearest, regrid , RMOC, remap

''' We have 9 runs:
    3,10,30,100,300,1000,3000, 10000 and Closed '''
if len(sys.argv) < 1:
   print '''\n    OOPSIE                                                                             
            usage:  give me a year to calculate over \n'''
   sys.exit(1)

Year = sys.argv[1]
tau=[3,10,30,100,300,1000,3000,10000,'Closed']

#Check what files exist
''' Let me know what files are missing and only bother
    if at least 3 exists'''
check=0
runs=[]
for i in range(len(tau)):
    flist=str(tau[i])+'daynokpp/'+Year+'all.nc'
    if not os.path.exists(flist):
       print 'WARNING: '+flist+' does not exist! (skipping this tau...)'
       check+=0
    else:
       check+=1
       runs.append(i)
if check<3:
   print "\n Less than 3 points to be plotted - no point. I'm stopping this..."
   sys.exit(1)

##~~~~~~~~~~~~~~~~Get the data~~~~~~~~~~~~~~~~~~~~~~###
Runs=np.array(runs)
#Grab a grid file to get model X,Y,Z arrays
file2=netcdf.netcdf_file(str(tau[Runs[0]])+'daynokpp/grid.nc','r')
Zp1=file2.variables['Zp1']
Zp=Zp1[:]*1
Z=file2.variables['Z']
Z=Z[:]
Y=file2.variables['Yp1']
Y=Y[:]*1
y=Y/1000
X=file2.variables['X']
X=X[:]*1
Yc=file2.variables['Y']
Yc=Yc[:]*1
yc=Yc/1000
Tau=[]
#General MOC
Psi=[]
#Take residuals 200k from Northern boundary
Psires1=[] #Upper 500m (AAIW)
Psires2=[] #500m - 1550m (NADW)
Psires3=[] #Below 1500m (AABW)
#Look at weird surface intensification
Surf=[] # 1200km upper 500m
# Look at AABW reduction near Souther boundary
AABW=[] # 200km from South
# Eddy maxima at 1000km
Eddy=[]
for i in range(len(Runs)):
    fname=str(tau[Runs[i]])+'daynokpp/'+Year+'all.nc'
    if tau[Runs[i]]=='Closed':
       Tau.append(100000)
    else:  
       Tau.append(tau[Runs[i]])
    #MOC
    psi,T=MOC(fname,Y,Z,Zp)
    d=np.nonzero(yc==find_nearest(yc,1200))[0][0]
    psimax=np.nanmax(psi[:,np.nonzero(y==find_nearest(y,1200))[0][0]]) 
    Psi.append(psimax)
    Rho = np.genfromtxt(str(tau[Runs[i]])+'daynokpp/Temp', delimiter = ',')
    Psimap,Psied,Zexp=remap(fname,Rho,X,Yc,Z,Y,Zp)
    #RMOC (Remapped)
    Q_levs = (np.arange(-10,10)+0.5)
    Psi_levs = Q_levs / 10
    Q_ticks = np.arange(-8,9,2.)
    Psi_ticks = Q_ticks / 10
    cf=plt.contourf(Yc,Zexp,Psimap,Psi_levs,cmap=plt.cm.seismic)
    cbar = plt.colorbar(cf, ticks=Psi_ticks, shrink=0.8)
    plt.clim(-1,1)
    plt.title("Remapped Stream Function ")
    plt.xlabel('Distance (km)')
    plt.ylabel('Depth (m)')
    cbar.ax.set_ylabel('$\psi_{res}$ (Sv)')
    plt.savefig("Psiremapped"+str(i)+".png")
    plt.clf()
    a=np.nonzero((Zexp)==find_nearest((Zexp),-500))[0][0]
    b=np.nonzero(yc==find_nearest(yc,yc[-1]-300))[0][0]
    c= np.nonzero(Zexp==find_nearest(Zexp,-1500))[0][0]
    e=np.nonzero(yc==find_nearest(yc,200))[0][0]
    psires1=np.nanmax(abs(Psimap[0:a,b]))
    psires2=np.nanmax(Psimap[a:c,b])
    psires3=np.nanmax(abs(Psimap[c::,b]))
    surf=np.nanmax(Psimap[:,d])   
    aabw=np.nanmax(abs(Psimap[:,e]))
    Psires1.append(psires1)
    Psires2.append(psires2)
    Psires3.append(psires3)
    Surf.append(surf)
    AABW.append(aabw)
    #Eddy 
    f=np.nonzero(yc==find_nearest(yc,1000))[0][0]
    eddy=np.nanmin(Psied[:,f])
    Eddy.append(eddy)
###~~~~~~~~~~~~~~~~~~~Make Plots~~~~~~~~~~~~~~~~~~~~~~###
#1) MOC
Relax=np.array(Tau)
Eulerian=np.array(Psi)
ax=plt.scatter(Relax,Eulerian,c='r', marker='o')
plt.title("Maximum Eulerian Mean MOC at 1200km (Averaged over years"+Year+")")
plt.xlabel(r'Relaxation timescale ($\tau$) in days')
plt.ylabel('MOC strength (Sv)')
plt.xscale('log')
plt.ylim(0,np.max(Eulerian)+0.5)
#ax.xaxis.set_ticks( Relax )
x=( os.path.expanduser('~')+"/Figures/")
plt.savefig(x+'/MOCALL'+Year+'.png')
plt.clf()
#2) RMOC MAIN
Psires1=np.array(Psires1)
Psires2=np.array(Psires2)
Psires3=np.array(Psires3)
fig1=plt.scatter(Relax,abs(Psires1),c='b', marker='o',label='AAIW')
fig2=plt.scatter(Relax,Psires2,c='r', marker='+',label='NADW')
fig3=plt.scatter(Relax,abs(Psires3),c='k', marker='^',label='AABW')
plt.title("Maximum RMOC at Northern edge of domain (Averaged over years"+Year+")")
plt.xlabel(r'Relaxation timescale ($\tau$) in days')
plt.ylabel('Absolute RMOC strength (Sv)')
plt.xscale('log')
plt.ylim(0,np.max(Psires2)+0.5)
plt.legend(loc='upper center', shadow=True)
x=( os.path.expanduser('~')+"/Figures/")
plt.savefig(x+'/RMOCALL'+Year+'.png')
plt.clf()
#3) RMOC Surface intensification
Surf=np.array(Surf)
plt.scatter(Relax,Surf,c='r', marker='o')
plt.title("Maximum RMOC at 1200km at 500m depth (Averaged over years"+Year+")")
plt.xlabel(r'Relaxation timescale ($\tau$) in days')
plt.ylabel('RMOC strength (Sv)')
plt.xscale('log')
plt.ylim(0,np.max(Surf)+0.1)
x=( os.path.expanduser('~')+"/Figures/")
plt.savefig(x+'/Surfaceintense'+Year+'.png')
plt.clf()
#4) RMOC AABW Reduction
AABW=np.array(AABW)
plt.scatter(Relax,abs(AABW),c='r', marker='o')
plt.title("Maximum AABW cell strength (Averaged over years "+Year+")")
plt.xlabel(r'Relaxation timescale ($\tau$) in days')
plt.ylabel('Absolute RMOC strength (Sv)')
plt.xscale('log')
plt.ylim(0,np.max(abs(AABW))+0.1)
x=( os.path.expanduser('~')+"/Figures/")
plt.savefig(x+'/AABWALL'+Year+'.png')
plt.clf()
#5) Eddy MOC 
Eddy=np.array(Eddy)
plt.scatter(Relax,abs(Eddy),c='r', marker='o')
plt.title("Maximum Eddy overturning  cell strength (1000km) (Averaged over years "+Year+")")
plt.xlabel(r'Relaxation timescale ($\tau$) in days')
plt.ylabel('Absolute Eddy MOC strength (Sv)')
plt.xscale('log')
plt.ylim(0,np.max(abs(Eddy))+0.5)
x=( os.path.expanduser('~')+"/Figures/")
plt.savefig(x+'/EddyALL'+Year+'.png')
plt.clf()
#6) ALL on ONE!!
plt.scatter(Relax,abs(Eddy),c='b', marker='o',label='Eddy')
plt.scatter(Relax,abs(AABW),c='b', marker='s',label='AABW(S)')
plt.scatter(Relax,abs(MOC),c='r', marker='x',label='MOC')
plt.scatter(Relax,abs(Sruf),c='r', marker='D',label='Surface MOC')
plt.scatter(Relax,abs(Psires1),c='b', marker='v',label='AAIW')
plt.scatter(Relax,Psires2,c='r', marker='+',label='NADW')
plt.scatter(Relax,abs(Psires3),c='b', marker='^',label='AABW')
plt.title("Max Stream Functions illustrating various cells (see reference figure)")
plt.xlabel(r'Relaxation timescale ($\tau$) in days')
plt.ylabel('RMOC strength (Sv)')
plt.xscale('log')
lgd=legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
#plt.ylim(0,np.max(Surf)+0.1)
x=( os.path.expanduser('~')+"/Figures/")
plt.savefig(x+'/AllStreamFunctions'+Year+'.png')
