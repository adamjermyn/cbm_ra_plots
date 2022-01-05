import mesa_reader as mr
import numpy as np
import matplotlib.gridspec as gridspec
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
prefix = '/Users/ajermyn/Dropbox/Active_Projects/CBM_trends/output/runs/'
DIR = prefix + 'main_Z_time_2021_12_14_13_00_27_sha_db89' + '/runs/' # The directory where you unpacked the data
FIGURES='./figures/' # Place to save plots
colors = [np.array(( 27,158,119, 255))/255,np.array((217, 95,  2, 255))/255,np.array((117,112,179, 255))/255]

def fit(m):
	a = 5.33591
	b = -0.823889
	c = 0.316251
	d = 0.152144
	e = -0.710093
	return np.sqrt(m-1.1)*np.tanh(m-1.1)*(a + b*m**2 + c*m**3+d*m**4) / (e + m**4)

def read_models(location,lis, z_getter, name, label):
	fig = plt.figure(figsize=(7,5))
	ax = plt.subplot(111)

	numcols, numrows = 200,200

	plt.gca().invert_xaxis()    

	x = []
	y = []
	zz = []
	zm = []
	zt = []

	for j in lis:
	  h=mr.MesaData(location+str(j)+'/LOGS/history.data')
	  model = h.model_number 
	  logl = h.log_L
	  logg = h.log_g 
	  loglh = h.log_LH 
	  center_h1 = h.center_h1 
	  logt= h.log_Teff 
	

	  ncz=h.subsurface_convective_regions 
	  tams=find_tams(center_h1,model)
	  zams=find_h(0.001,center_h1,model)
	  mams=find_mams(center_h1,model)
		
	  # Create Lists  
	  zz.append(z_getter(h)[zams])
	  zm.append(z_getter(h)[mams])
	  zt.append(z_getter(h)[tams])

	fig = plt.figure(figsize=(6,5))
	spec = gridspec.GridSpec(ncols=2,nrows=1, figure=fig, width_ratios=(1,1), wspace=0, hspace=0)
	axes = list(fig.add_subplot(spec[j]) for j in range(2))

	fitted = fit(np.array(lis))

	axes[0].plot(lis, zz, color=colors[0], label='ZAMS ($X_c=0.72$)')
	axes[0].plot(lis, zm, color=colors[1], label='mid-MS ($X_c = 0.36$)')
	axes[0].plot(lis, zt, color=colors[2], label='TAMS ($X_c = 0$)')
	axes[0].plot(lis, fitted, color=colors[1], linestyle=':', label='Fitted mid-MS')
	axes[0].set_xlabel(r'$M/M_\odot$')
	axes[0].set_ylabel(label)
	axes[0].set_xlim([1,3])
	axes[1].plot(lis, zz, color=colors[0], label='ZAMS ($X_c=0.72$)')
	axes[1].plot(lis, zm, color=colors[1], label='mid-MS ($X_c = 0.36$)')
	axes[1].plot(lis, zt, color=colors[2], label='TAMS ($X_c = 0$)')
	axes[1].plot(lis, fitted, color=colors[1], linestyle=':', label='Fitted mid-MS')
	axes[1].legend()
	axes[1].set_xlabel(r'$M/M_\odot$')
	axes[1].set_xlim([3,60])

	print(zm)

	# Remove inner tick labels
	axes[1].set_yticklabels([])

	plt.savefig(FIGURES+name,bbox_inches='tight')

read_models(DIR,mods,dr_core_div_h_getter, 'fov_mass.pdf', r'$\alpha_{\rm ov}$')
