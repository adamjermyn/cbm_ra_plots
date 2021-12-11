import mesa_reader as mr
import numpy as np
from pylab import *
from lowPr_highTa_onset import Ra_crit
from math import log10, pi

from scipy import interpolate
from scipy.interpolate import griddata

from numpy import loadtxt
import pandas as pd 
 
from constants import *
from plot_settings import *
from functions import *
from getters import *

mods = [1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10,10.5,11,11.5,12,14,16,18,20,21,22,23,24,25,26,27,28,29,30,32,35,37,40,42,45,47,50,52,55,57,60]
hrdlines = [1.5,2.0,2.4,3.0,4.0,5.0,6.0,7.0,9.0,12,20,30,40,60]

# Data Location
FIGURES='./figures/' # Place to save plots
prefix = '/Users/ajermyn/Dropbox/Active_Projects/CBM_trends/output/runs/'

# STRINGS
logteff=r'$\log_{10}\, T_{\rm eff}$/K'
logell=r'$\log_{10}\, L$/L$_\odot$'

Ta = 0

def tri_area(xs,ys):
  arr = np.ones((3,3))
  arr[0] = xs
  arr[1] = ys
  area = 0.5 * np.linalg.det(arr)
  return area

def read_models(location,lis, z_getter, Pr_getter, name, ann, bar_label, cmap, extra_label):
    fig = plt.figure(figsize=(7,5))
    ax = plt.subplot(111)

    numcols, numrows = 200,200

    plt.gca().invert_xaxis()    

    x = []
    y = []
    z = []

    for j in lis:
      h=mr.MesaData(location+str(j)+'/LOGS/history.data')
      model = h.model_number 
      logl = h.log_L
      logg = h.log_g 
      loglh = h.log_LH 
      center_h1 = h.center_h1 
      logt= h.log_Teff 
    

      ncz=h.subsurface_convective_regions 
      zams=find_zams(logl,loglh,model)
      tams=find_tams(center_h1,model)
      zams=find_h(0.001,center_h1,model)
        
      # Create Lists  
      x.append(logt[zams:])
      y.append(logl[zams:])
      z.append(z_getter(h)[zams:])
      Pr = 10**Pr_getter(h)[zams:]

      for i in range(len(Pr)):
        z[-1][i] -= np.log10(Ra_crit(Ta, Pr[i]))

    x=array(list(flatten(x)))
    y=array(list(flatten(y)))
    z=array(list(flatten(z)))

    xi = np.linspace(x.min(), x.max(), numcols)
    yi = np.linspace(y.min(), y.max(), numrows)
    
    triang = tri.Triangulation(x,y)
    areas = np.array(list(tri_area(x[triang.triangles[q]], y[triang.triangles[q]]) for q in range(len(triang.triangles))))
    triang.set_mask((areas > 0.005))

    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    ax.contour(xi, yi, zi, 9, levels=[-8,-6,-4,-2,0,2,4,6,8], extend='both', colors='k')    
    cntr1 = ax.contourf(xi, yi, zi, 9, levels=[-8,-6,-4,-2,0,2,4,6,8], extend='both', cmap=cmap)
    cbar = fig.colorbar(cntr1, ax=ax)
    cbar.ax.set_ylabel(bar_label)
    ls = list(cbar.ax.get_yticks())
    cbar.ax.set_yticklabels(['{:.2f}'.format(x) for x in ls])

    dt = 0.07
    dl = -0.09
    
    for j in hrdlines:
      h=mr.MesaData(DIR+str(j)+'/LOGS/history.data')
      model = h.model_number 
      logl = h.log_L
      logg = h.log_g 
      loglh=h.log_LH 
      center_h1 = h.center_h1 
      logt= h.log_Teff 
      zams=find_zams(logl,loglh,model)
      tams=find_tams(center_h1,model)
      zams=find_h(0.001,center_h1,model)
      ax.plot(logt[zams:],logl[zams:],c='gray',alpha=1.0)

      if j == 12:
        # Fixes up the placement so 11Msun and 12Msun don't run into each other.
        ax.text(logt[zams]+dt,logl[zams]+dl+0.1,str(j)+r'$M_\odot$',ha='center',fontsize=14) # ,verticalalignment='center',rotation='vertical',  
      else:
        ax.text(logt[zams]+dt,logl[zams]+dl,str(j)+r'$M_\odot$',ha='center',fontsize=14) # ,verticalalignment='center',rotation='vertical',  
                
    ax.set_xlabel(logteff)
    ax.set_ylabel(logell)
    ax.text(3.8,5.5,ann,ha='center',fontsize=18)
    ax.set_xlim([4.83,3.7])
    ax.set_ylim([0,6.1])   
    plt.savefig(FIGURES+name,bbox_inches='tight')


DIR = prefix + 'main_Z_time_2021_12_09_17_42_58_sha_dab3' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,Ra_HI_getter, Pr_HI_getter, 'HI_Ra_Z_Z_MW.pdf', 'HI', r'$\log\mathrm{Ra}/\mathrm{Ra}_{\rm crit}$', cmap='PiYG', extra_label=None)
read_models(DIR,mods,Ra_HeI_getter, Pr_HeI_getter, 'HeI_Ra_Z_Z_MW.pdf', 'HeI', r'$\log\mathrm{Ra}/\mathrm{Ra}_{\rm crit}$', cmap='PiYG', extra_label=None)
read_models(DIR,mods,Ra_HeII_getter, Pr_HeII_getter, 'HeII_Ra_Z_Z_MW.pdf', 'HeII', r'$\log\mathrm{Ra}/\mathrm{Ra}_{\rm crit}$', cmap='PiYG', extra_label=None)
read_models(DIR,mods,Ra_FeCZ_getter, Pr_FeCZ_getter, 'FeCZ_Ra_Z_Z_MW.pdf', 'FeCZ', r'$\log\mathrm{Ra}/\mathrm{Ra}_{\rm crit}$', cmap='PiYG', extra_label=None)
