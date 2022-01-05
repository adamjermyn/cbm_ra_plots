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
hrdlines = [7.0,8.0,9.0,11,16,20,25,30,40]

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

def read_models(location,lis, period, z_getters, Pr_getters, name, ann, extra_label):
    fig = plt.figure(figsize=(7,5))
    ax = plt.subplot(111)

    numcols, numrows = 200,200

    plt.gca().invert_xaxis()    

    omega = 2 * np.pi / (24 * 3600 * period) # Convert period in days to Omega in rad/s

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

      # Correct for Taylor-number scaling
      viscous_times = viscous_times_getter(h)[zams:,:]
      Ta = 4 * (viscous_times * omega)**2

      for z_getter,Pr_getter,k in zip(*(z_getters,Pr_getters,range(4))):
        z[-1].append(z_getter(h)[zams:])

        Pr = 10**Pr_getter(h)[zams:]

        for i in range(len(Pr)):
          z[-1][-1][i] -= np.log10(Ra_crit(Ta[i,k], Pr[i]))

    for i in range(len(z)):
      z[i] = np.array(z[i])
      z[i][np.isnan(z[i])] = -np.inf
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

    cmap3 = CustomCmap([0.02, 0.75, 1.00], [0.02, 0.75, 1.00]) # from white to +/- 5,192,255
    ax.pcolormesh(xi, yi, zi, cmap=cmap3, alpha=0.4, edgecolors='face')
    
    dt = 0.04
    dl = -0.06
    
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
    ax.text(4.74,2.35,ann,ha='left',fontsize=18)
    ax.set_xlim([4.75,4.3])
    ax.set_ylim([2.3,4])   
    plt.savefig(FIGURES+name,bbox_inches='tight')

Ra_getters = [Ra_HI_getter,Ra_HeI_getter,Ra_HeII_getter,Ra_FeCZ_getter]
Pr_getters = [Pr_HI_getter,Pr_HeI_getter,Pr_HeII_getter,Pr_FeCZ_getter]

DIR = prefix + 'main_Z_time_2021_12_14_12_59_57_sha_3d07' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,1000,Ra_getters,Pr_getters, 'spec_no_cz_Z_Z_SMC_P_1000.pdf', r'$Z=Z_{\rm SMC}=0.002$', extra_label=None)

DIR = prefix + 'main_Z_time_2021_12_14_13_00_09_sha_5afc' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,1000,Ra_getters,Pr_getters, 'spec_no_cz_Z_Z_LMC_P_1000.pdf', r'$Z=Z_{\rm LMC}=0.006$', extra_label=None)

DIR = prefix + 'main_Z_time_2021_12_14_13_00_18_sha_c740' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,1000,Ra_getters,Pr_getters, 'spec_no_cz_Z_0.01_P_1000.pdf', r'$Z=0.01$', extra_label=r'$P=1000\,\mathrm{d}$')

DIR = prefix + 'main_Z_time_2021_12_14_13_00_27_sha_db89' + '/runs/' # The directory where you unpacked the data
read_models(DIR,mods,1000,Ra_getters,Pr_getters, 'spec_no_cz_Z_Z_MW_P_1000.pdf', r'$Z=Z_{\rm MW}=0.014$', extra_label=None)
