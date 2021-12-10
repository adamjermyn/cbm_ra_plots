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

mods = [7.0,7.5,8.0,8.5,9.0,9.5,10,10.5,11,11.5,12,14,16,18,20,21,22,23,24,25,26,27,28,29,30,32,35,37,40]
hrdlines = [8.0,9.0,12,16,20,25,30,40]

# Data Location
FIGURES='./figures/' # Place to save plots
prefix = '/Users/ajermyn/Dropbox/Active_Projects/CBM_trends/output/runs/'

# STRINGS
logteff=r'$\log_{10}\, T_{\rm eff}$/K'
logell=r'$\log_{10}\, \mathscr{L}/\mathscr{L}_\odot$'

def tri_area(xs,ys):
  arr = np.ones((3,3))
  arr[0] = xs
  arr[1] = ys
  area = 0.5 * np.linalg.det(arr)
  return area

def read_models(location,lis, z_getters, name, ann, extra_label):
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
      y.append(logl[zams:]-np.log10(j))        
      z.append([])
      for z_getter,Pr_getter in zip(*(z_getters,Pr_getters)):
        z[-1].append(z_getter(h)[zams:])

        Pr = Pr_getter(h)[zams:]

        for i in range(len(Pr)):
          z[-1][-1][i] -= np.log10(Ra_crit(Ta, Pr[i]))

    for i in range(len(z)):
      z[i] = np.array(z[i])
      z[i] = np.amax(z[i], axis=0)
      z[i][z[i] > 0] = np.nan
      z[i][z[i] < 0] = 1

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

    ax.pcolormesh(xi, yi, zi, edgecolors='face')
    
    dt = 0.03
    dl = -0.05
    
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
      ax.plot(logt[zams:],logl[zams:]-np.log10(j),c='gray',alpha=1.0)

      if j == 12:
        # Fixes up the placement so 11Msun and 12Msun don't run into each other.
        ax.text(logt[zams]+dt,logl[zams]-np.log10(j)+dl+0.1,str(j)+r'$M_\odot$',ha='center',fontsize=14) # ,verticalalignment='center',rotation='vertical',  
      else:
        ax.text(logt[zams]+dt,logl[zams]-np.log10(j)+dl,str(j)+r'$M_\odot$',ha='center',fontsize=14) # ,verticalalignment='center',rotation='vertical',  
                
    ax.set_xlabel(logteff)
    ax.set_ylabel(logell)
    ax.text(4.74,2.45,ann,ha='left',fontsize=18)
    ax.set_xlim([4.75,4.3])
    ax.set_ylim([2.4,4])   
    plt.savefig(FIGURES+name,bbox_inches='tight')

Ra_getters = [Ra_HI_getter,Ra_HeI_getter,Ra_HeII_getter,Ra_FeCZ_getter]
Pr_getters = [Pr_HI_getter,Pr_HeI_getter,Pr_HeII_getter,Pr_FeCZ_getter]

DIR = prefix + 'main_Z_time_2021_12_09_17_42_33_sha_db8c' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,Ra_getters,Pr_getters, 'spec_no_cz_Z_Z_SMC.pdf', r'$Z=Z_{\rm SMC}=0.002$', extra_label=None)

DIR = prefix + 'main_Z_time_2021_12_09_17_42_41_sha_0199' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,Ra_getters,Pr_getters, 'spec_no_cz_Z_Z_LMC.pdf', r'$Z=Z_{\rm LMC}=0.006$', extra_label=None)

DIR = prefix + 'main_Z_time_2021_12_09_17_42_49_sha_1dcf' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,Ra_getters,Pr_getters, 'spec_no_cz_Z_0.01.pdf', r'$Z=0.01$', extra_label=None)

DIR = prefix + 'main_Z_time_2021_12_09_17_42_58_sha_dab3' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,Ra_getters,Pr_getters, 'spec_no_cz_Z_Z_MW.pdf', r'$Z=Z_{\rm MW}=0.014$', extra_label=None)
