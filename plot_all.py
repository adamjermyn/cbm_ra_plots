import mesa_reader as mr
import numpy as np
from pylab import *
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
DIR = prefix + 'main_Z_time_2021_12_08_13_30_07_sha_2c98' + '/runs/' # The directory where you unpacked the data

# STRINGS
logteff=r'$\log_{10}\, T_{\rm eff}$/K'
logell=r'$\log_{10}\, L$/L$_\odot$'

# Use ell_sun to make Spectroscopic HRD, if needed
ell_sun=(5777)**4.0/(274*100)  

def tri_area(xs,ys):
  arr = np.ones((3,3))
  arr[0] = xs
  arr[1] = ys
  area = 0.5 * np.linalg.det(arr)
  return area

def read_models(location,lis, z_getter, name, bar_label, cmap, extra_label):
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

    x=array(list(flatten(x)))
    y=array(list(flatten(y)))
    z=array(list(flatten(z)))

    if name == 'dlnm_core_one_panel.pdf':
      z[z > 0.6] = np.nan
    if name == 'dlnr_core_one_panel.pdf':
      z[z > 0.4] = np.nan
    if name == 'Ra.pdf':
      z[z < 0] = np.nan

    xi = np.linspace(x.min(), x.max(), numcols)
    yi = np.linspace(y.min(), y.max(), numrows)
    
    triang = tri.Triangulation(x,y)
    areas = np.array(list(tri_area(x[triang.triangles[q]], y[triang.triangles[q]]) for q in range(len(triang.triangles))))
    triang.set_mask((areas > 0.005))

    interpolator = tri.LinearTriInterpolator(triang, z)
    Xi, Yi = np.meshgrid(xi, yi)
    zi = interpolator(Xi, Yi)

    ax.contour(xi, yi, zi, 8, colors='k')    
    cntr1 = ax.contourf(xi, yi, zi, 8, cmap=cmap)
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
      ell = (10**logt)**4.0/(10**logg)
      ell=np.log10(ell/ell_sun)  
      ax.plot(logt[zams:],logl[zams:],c='gray',alpha=1.0)

      if j == 12:
        # Fixes up the placement so 11Msun and 12Msun don't run into each other.
        ax.text(logt[zams]+dt,logl[zams]+dl+0.1,str(j)+r'$M_\odot$',ha='center',fontsize=14) # ,verticalalignment='center',rotation='vertical',  
      else:
        ax.text(logt[zams]+dt,logl[zams]+dl,str(j)+r'$M_\odot$',ha='center',fontsize=14) # ,verticalalignment='center',rotation='vertical',  
                
    ax.set_xlabel(logteff)
    ax.set_ylabel(logell)
    
    ax.set_xlim([4.83,3.7])
    ax.set_ylim([0,6.1])   
    plt.savefig(FIGURES+name,bbox_inches='tight')


read_models(DIR,mods,original_m_core_over_m,'original_m_core_over_m.pdf',r'$M_{\mathrm{core}}/M_\star$',cmap='RdPu',extra_label=None)
read_models(DIR,mods,m_core_over_m,'m_core_over_m.pdf',r'$M_{\mathrm{core}}/M_\star$',cmap='RdPu',extra_label=None)
read_models(DIR,mods,dr_core_div_h_getter, 'dr_core_div_h_core_one_panel.pdf', r'$\alpha_{\rm ov} = \Delta R_{\mathrm{core}}/h$', cmap='YlOrRd', extra_label=None)
read_models(DIR,mods,dlnm_core_getter, 'dlnm_core_one_panel.pdf', r'$\Delta M_{\mathrm{core}}/M_{\mathrm{core}}$', cmap='BuPu', extra_label=None)
read_models(DIR,mods,dlnr_core_getter, 'dlnr_core_one_panel.pdf', r'$\Delta R_{\mathrm{core}}/R_{\mathrm{core}}$', cmap='YlGnBu', extra_label=None)

exit()


read_models(DIR,mods,Re_HI_getter, 'HI_Re.pdf', r'HI $\log\mathrm{Re}$', cmap='YlGnBu', extra_label=None)
read_models(DIR,mods,Pr_HI_getter, 'HI_Pr.pdf', r'HI $\log\mathrm{Pr}$', cmap='BuPu', extra_label=None)

read_models(DIR,mods,Re_HeI_getter, 'HeI_Re.pdf', r'HeI $\log\mathrm{Re}$', cmap='YlGnBu', extra_label=None)
read_models(DIR,mods,Pr_HeI_getter, 'HeI_Pr.pdf', r'HeI $\log\mathrm{Pr}$', cmap='BuPu', extra_label=None)

read_models(DIR,mods,Re_HeII_getter, 'HeII_Re.pdf', r'HeII $\log\mathrm{Re}$', cmap='YlGnBu', extra_label=None)
read_models(DIR,mods,Pr_HeII_getter, 'HeII_Pr.pdf', r'HeII $\log\mathrm{Pr}$', cmap='BuPu', extra_label=None)

read_models(DIR,mods,Re_FeCZ_getter, 'FeCZ_Re.pdf', r'FeCZ $\log\mathrm{Re}$', cmap='YlGnBu', extra_label=None)
read_models(DIR,mods,Pr_FeCZ_getter, 'FeCZ_Pr.pdf', r'FeCZ $\log\mathrm{Pr}$', cmap='BuPu', extra_label=None)




read_models(DIR,mods,Ra_getter, 'Ra.pdf', r'$\log\mathrm{Ra}$', cmap='YlOrRd', extra_label=None)
read_models(DIR,mods,Re_getter, 'Re.pdf', r'$\log\mathrm{Re}$', cmap='YlGnBu', extra_label=None)
read_models(DIR,mods,Pr_getter, 'Pr.pdf', r'$\log\mathrm{Pr}$', cmap='BuPu', extra_label=None)


